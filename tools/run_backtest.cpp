// Runs the Phase 15 historical delta-hedging backtest (backtest.hpp) for
// all three pricing models -- Black-Scholes, Heston, and Bates -- over
// the full BTC price history in data/btc_price_history.csv, and emits a
// JSON summary to stdout. The same "never hand-type real numbers into
// docs" convention as gen_stats.cpp and calibrate_heston.cpp: this
// program is the source of truth for every backtest number that ends up
// in README.md/docs/numerics.md/docs/phase-log.md.
//
// Heston and Bates structural parameters below are NOT fit inside this
// program -- they're the same real, already-calibrated values this
// project's other tools produce:
//   - Heston kappa/theta/xi/rho: tools/calibrate_heston.cpp's output,
//     fit to the real committed Deribit snapshot (data/btc_chain_snapshot.json).
//   - Bates jump_intensity/jump_mean/jump_vol: estimated directly from
//     this same price history via a simple 4-standard-deviation
//     threshold on daily log returns (see docs/numerics.md Phase 15 for
//     the method and its honest limitations) -- computed once with a
//     throwaway script, not re-derived here, and hardcoded the same way
//     the Heston values are.
// Both are held fixed across the entire ~16-year backtest; only each
// window's v0 (instantaneous variance) is re-estimated per window, from
// trailing realized volatility. See backtest.hpp's header comment for
// why.

#include "deriv-engine/backtest.hpp"
#include "deriv-engine/price_history.hpp"
#include "deriv-engine/black_scholes.hpp"
#include "deriv-engine/greeks.hpp"
#include "deriv-engine/heston.hpp"
#include "deriv-engine/bates.hpp"
#include <cstdio>
#include <cmath>

using namespace deriv;

namespace {

// Central finite-difference delta: the only option available for
// Heston/Bates, which have no analytical Greeks in this engine (see
// backtest.hpp's header comment for why that's a deliberate scope
// decision, not an oversight). h is a relative bump (0.1% of spot) --
// small enough to approximate the true derivative closely, large enough
// to stay well clear of the COS method's own numerical noise floor at
// double precision.
double finite_diff_delta(double spot, std::function<double(double)> price_at_spot) {
    double h = spot * 0.001;
    return (price_at_spot(spot + h) - price_at_spot(spot - h)) / (2.0 * h);
}

ModelFactory make_black_scholes_factory() {
    return [](double strike, double vol) -> PriceDeltaFn {
        return [strike, vol](double spot, double tau) -> std::pair<double, double> {
            EuropeanOption option{strike, tau, OptionType::Call};
            MarketData market{spot, 0.0, 0.0, vol};
            double price = black_scholes_price(option, market);
            Greeks g = black_scholes_greeks(option, market);
            return {price, g.delta};
        };
    };
}

ModelFactory make_heston_factory(double kappa, double theta_structural, double xi, double rho) {
    // theta_structural is the real calibrated long-run variance; each
    // window's v0 is that window's own trailing realized variance
    // instead (see backtest.hpp) -- theta stays fixed as the structural
    // "where variance reverts to" parameter regardless.
    return [kappa, theta_structural, xi, rho](double strike, double window_vol) -> PriceDeltaFn {
        double v0 = window_vol * window_vol;
        HestonParams params{v0, kappa, theta_structural, xi, rho};
        return [strike, params](double spot, double tau) -> std::pair<double, double> {
            EuropeanOption option{strike, tau, OptionType::Call};
            auto price_at = [&](double s) {
                MarketData m{s, 0.0, 0.0, 0.0};
                return heston_price(option, m, params);
            };
            double price = price_at(spot);
            double delta = finite_diff_delta(spot, price_at);
            return {price, delta};
        };
    };
}

ModelFactory make_bates_factory(double kappa, double theta_structural, double xi, double rho, double jump_intensity,
                                 double jump_mean, double jump_vol) {
    return [=](double strike, double window_vol) -> PriceDeltaFn {
        double v0 = window_vol * window_vol;
        BatesParams params{{v0, kappa, theta_structural, xi, rho}, jump_intensity, jump_mean, jump_vol};
        return [strike, params](double spot, double tau) -> std::pair<double, double> {
            EuropeanOption option{strike, tau, OptionType::Call};
            auto price_at = [&](double s) {
                MarketData m{s, 0.0, 0.0, 0.0};
                return bates_price(option, m, params);
            };
            double price = price_at(spot);
            double delta = finite_diff_delta(spot, price_at);
            return {price, delta};
        };
    };
}

void print_summary(const BacktestSummary& s) {
    std::printf("  {\n");
    std::printf("    \"model\": \"%s\",\n", s.model_name.c_str());
    std::printf("    \"num_windows\": %d,\n", s.num_windows);
    std::printf("    \"mean_pnl_bps\": %.4f,\n", s.mean_pnl_bps);
    std::printf("    \"stdev_pnl_bps\": %.4f\n", s.stdev_pnl_bps);
    std::printf("  }");
}

}  // namespace

int main(int argc, char** argv) {
    auto history = load_price_history("data/btc_price_history.csv");

    // Optional --since=YYYY-MM-DD to trim the leading (extremely thin,
    // extremely volatile) years of BTC's trading history -- a sensitivity
    // check, not the headline number. See docs/numerics.md Phase 15 for
    // why: 2010-2013 prices move on razor-thin volume, and Deribit
    // itself (the source of every OTHER real-market number in this
    // project) didn't exist before 2016, so "was this ever a market
    // options could realistically trade in" is a fair question to ask
    // separately from "does the full numerical history support this
    // conclusion."
    std::string since;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        const std::string prefix = "--since=";
        if (arg.rfind(prefix, 0) == 0) since = arg.substr(prefix.size());
    }
    if (!since.empty()) {
        std::vector<PricePoint> filtered;
        for (const auto& p : history) {
            if (p.date >= since) filtered.push_back(p);
        }
        history = filtered;
    }

    std::printf("{\n");
    std::printf("  \"history_points\": %d,\n", static_cast<int>(history.size()));
    std::printf("  \"date_range\": [\"%s\", \"%s\"],\n", history.front().date.c_str(), history.back().date.c_str());
    std::printf("  \"results\": [\n");

    // Heston structural params: tools/calibrate_heston.cpp's real output
    // against data/btc_chain_snapshot.json (99-day expiry, 66 contracts,
    // converged, 0.43pp RMSE) -- kappa=1.829505, theta=0.134505,
    // xi=1.280247, rho=-0.194560.
    const double kappa = 1.829505, theta = 0.134505, xi = 1.280247, rho = -0.194560;

    // Bates jump params: estimated from this same price history via the
    // 4-sigma threshold method (docs/numerics.md Phase 15).
    const double jump_intensity = 3.4626461121463663;
    const double jump_mean = -0.000773771446447099;
    const double jump_vol = 0.27331845903521157;

    auto bs = run_backtest(history, "Black-Scholes", make_black_scholes_factory());
    print_summary(bs);
    std::printf(",\n");

    auto heston = run_backtest(history, "Heston", make_heston_factory(kappa, theta, xi, rho));
    print_summary(heston);
    std::printf(",\n");

    auto bates = run_backtest(history, "Bates",
                               make_bates_factory(kappa, theta, xi, rho, jump_intensity, jump_mean, jump_vol));
    print_summary(bates);
    std::printf("\n");

    std::printf("  ]\n");
    std::printf("}\n");

    return 0;
}
