#include "deriv-engine/var.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

namespace deriv {

namespace {

// Empirical quantile via the nearest-rank method on a distribution
// that's already sorted ascending: the value such that a fraction `p` of
// the sample is at or below it. Not interpolated between ranks -- simple
// and standard for VaR reporting, and the sample sizes used here
// (thousands of Monte Carlo paths, or the full historical return
// series) make the difference from a smoother interpolation negligible.
double empirical_quantile(const std::vector<double>& sorted_ascending, double p) {
    if (sorted_ascending.empty()) throw std::runtime_error("empirical_quantile: empty sample");
    int n = static_cast<int>(sorted_ascending.size());
    int idx = static_cast<int>(std::floor(p * n));
    idx = std::clamp(idx, 0, n - 1);
    return sorted_ascending[idx];
}

VaRResult summarize_pnl(const std::string& method_name, std::vector<double> pnl) {
    std::sort(pnl.begin(), pnl.end());

    double sum = 0.0;
    for (double x : pnl) sum += x;
    double mean = sum / pnl.size();

    // VaR is reported as a positive loss number: the 5th/1st percentile
    // of the (signed) P&L distribution is itself a large negative
    // number for a risky book, so negate it.
    double var_95 = -empirical_quantile(pnl, 0.05);
    double var_99 = -empirical_quantile(pnl, 0.01);

    return VaRResult{method_name, static_cast<int>(pnl.size()), mean, var_95, var_99};
}

// One path's terminal spot after `horizon_days` daily full-truncation
// Euler steps under the Bates SDE (Heston's variance process, exactly
// as in heston_mc.cpp/bates_mc.cpp, plus a compound Poisson jump added
// to the log-price increment each step -- pass jump_intensity=0 in
// `params` for a pure-Heston simulation). This is a THIRD independent
// implementation of the same stepping scheme heston_mc.cpp and
// bates_mc.cpp each already have their own copy of -- consistent with
// this project's existing precedent (those two files don't share their
// stepping logic with each other either) rather than a new pattern
// introduced here.
double simulate_bates_terminal_spot(double S0, double r, double q, const BatesParams& params, int horizon_days,
                                     std::mt19937& rng, std::normal_distribution<double>& normal) {
    const double dt = 1.0 / 365.0;
    const double sqrt_dt = std::sqrt(dt);
    const auto& hp = params.heston;
    const double sqrt_1_minus_rho_sq = std::sqrt(std::max(1.0 - hp.rho * hp.rho, 0.0));

    double k = std::exp(params.jump_mean + 0.5 * params.jump_vol * params.jump_vol) - 1.0;
    double r_compensated = r - params.jump_intensity * k;

    std::poisson_distribution<int> poisson(params.jump_intensity * dt);

    double v = hp.v0;
    double log_S = std::log(S0);

    for (int step = 0; step < horizon_days; ++step) {
        double v_plus = std::max(v, 0.0);
        double sqrt_v_plus = std::sqrt(v_plus);

        double z1 = normal(rng);
        double z2 = normal(rng);
        double z_var = hp.rho * z1 + sqrt_1_minus_rho_sq * z2;

        double jump_sum = 0.0;
        int num_jumps = poisson(rng);
        for (int j = 0; j < num_jumps; ++j) {
            jump_sum += params.jump_mean + params.jump_vol * normal(rng);
        }

        log_S += (r_compensated - q - 0.5 * v_plus) * dt + sqrt_v_plus * sqrt_dt * z1 + jump_sum;
        v = v + hp.kappa * (hp.theta - v_plus) * dt + hp.xi * sqrt_v_plus * sqrt_dt * z_var;
    }

    return std::exp(log_S);
}

}  // namespace

VaRResult historical_var(const Portfolio& portfolio, const MarketData& market,
                          const std::vector<PricePoint>& history, const PortfolioPricer& pricer,
                          int horizon_days) {
    if (static_cast<int>(history.size()) <= horizon_days) {
        throw std::runtime_error("historical_var: not enough history for this horizon");
    }

    double base_value = portfolio_value(portfolio, market, pricer);
    Portfolio decayed = decay_time(portfolio, horizon_days);

    std::vector<double> pnl;
    int n = static_cast<int>(history.size());
    pnl.reserve(n - horizon_days);

    for (int i = 0; i + horizon_days < n; ++i) {
        double log_return = std::log(history[i + horizon_days].price / history[i].price);
        double shocked_spot = market.spot * std::exp(log_return);

        MarketData shocked_market = market;
        shocked_market.spot = shocked_spot;

        double shocked_value = portfolio_value(decayed, shocked_market, pricer);
        pnl.push_back(shocked_value - base_value);
    }

    return summarize_pnl("historical", pnl);
}

VaRResult monte_carlo_var(const Portfolio& portfolio, const MarketData& market, const BatesParams& params,
                           const PortfolioPricer& pricer, int num_paths, int horizon_days, unsigned int seed) {
    double base_value = portfolio_value(portfolio, market, pricer);
    Portfolio decayed = decay_time(portfolio, horizon_days);

    std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
    std::normal_distribution<double> normal(0.0, 1.0);

    std::vector<double> pnl;
    pnl.reserve(num_paths);

    for (int i = 0; i < num_paths; ++i) {
        double shocked_spot =
            simulate_bates_terminal_spot(market.spot, market.risk_free_rate, market.dividend_yield, params,
                                          horizon_days, rng, normal);

        MarketData shocked_market = market;
        shocked_market.spot = shocked_spot;

        double shocked_value = portfolio_value(decayed, shocked_market, pricer);
        pnl.push_back(shocked_value - base_value);
    }

    return summarize_pnl("monte_carlo", pnl);
}

double delta_normal_var(const Portfolio& portfolio, const MarketData& market, const PortfolioPricer& pricer,
                         double annualized_vol, int horizon_days, double confidence) {
    double delta = portfolio_delta(portfolio, market, pricer);

    double z;
    if (std::abs(confidence - 0.95) < 1e-9) {
        z = 1.644854;
    } else if (std::abs(confidence - 0.99) < 1e-9) {
        z = 2.326348;
    } else {
        throw std::runtime_error("delta_normal_var: only 0.95 and 0.99 confidence are supported");
    }

    double horizon_years = horizon_days / 365.0;
    double move_size = market.spot * annualized_vol * std::sqrt(horizon_years);

    // Linear approximation: assumes the portfolio's P&L is delta * (spot
    // move), ignoring gamma/vega/every other nonlinearity entirely --
    // see var.hpp's header comment.
    return z * std::abs(delta) * move_size;
}

StressResult run_stress_test(const Portfolio& portfolio, const MarketData& market,
                              const std::vector<PricePoint>& history, const PortfolioPricer& pricer,
                              const StressScenario& scenario) {
    auto find_date = [&](const std::string& date) -> const PricePoint& {
        auto it = std::lower_bound(history.begin(), history.end(), date,
                                    [](const PricePoint& p, const std::string& d) { return p.date < d; });
        if (it == history.end() || it->date != date) {
            throw std::runtime_error("run_stress_test: date not found in history: " + date);
        }
        return *it;
    };

    const PricePoint& start = find_date(scenario.start_date);
    const PricePoint& end = find_date(scenario.end_date);

    double pct_change = end.price / start.price - 1.0;

    // Actual number of calendar days between the two dates, read off
    // their position in `history` (which is daily and gap-free -- see
    // price_history.hpp) rather than parsed from the date strings, to
    // avoid pulling in calendar-arithmetic code this project doesn't
    // otherwise need.
    int start_idx = static_cast<int>(&start - &history[0]);
    int end_idx = static_cast<int>(&end - &history[0]);
    int elapsed_days = end_idx - start_idx;

    double base_value = portfolio_value(portfolio, market, pricer);
    Portfolio decayed = decay_time(portfolio, elapsed_days);

    MarketData shocked_market = market;
    shocked_market.spot = market.spot * (1.0 + pct_change);

    double shocked_value = portfolio_value(decayed, shocked_market, pricer);

    return StressResult{scenario.label,       scenario.start_date, scenario.end_date, start.price,
                         end.price,            pct_change,          base_value,        shocked_value,
                         shocked_value - base_value};
}

std::vector<StressScenario> default_stress_scenarios() {
    return {
        {"2020 COVID crash", "2020-02-19", "2020-03-13"},
        {"2022 FTX collapse", "2022-11-06", "2022-11-21"},
        {"2021-2022 bear market (peak to trough)", "2021-11-10", "2022-11-21"},
    };
}

}  // namespace deriv
