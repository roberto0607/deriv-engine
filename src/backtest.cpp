#include "deriv-engine/backtest.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace deriv {

double trailing_realized_vol(const std::vector<PricePoint>& history, int start_idx, int num_days) {
    if (start_idx < 0 || num_days <= 0 || start_idx + num_days >= static_cast<int>(history.size())) {
        throw std::runtime_error("trailing_realized_vol: window out of range");
    }

    std::vector<double> log_returns;
    log_returns.reserve(num_days);
    for (int i = 0; i < num_days; ++i) {
        double s0 = history[start_idx + i].price;
        double s1 = history[start_idx + i + 1].price;
        log_returns.push_back(std::log(s1 / s0));
    }

    double mean = std::accumulate(log_returns.begin(), log_returns.end(), 0.0) / log_returns.size();
    double sq_diff_sum = 0.0;
    for (double r : log_returns) {
        double d = r - mean;
        sq_diff_sum += d * d;
    }
    // Sample (n-1) standard deviation -- num_days is typically 30, small
    // enough that the n vs n-1 distinction isn't negligible.
    double sample_variance = sq_diff_sum / (log_returns.size() - 1);
    double daily_vol = std::sqrt(sample_variance);

    // 365, not the ~252 trading-day convention equities use -- BTC
    // trades (and this history has a price point) every calendar day.
    return daily_vol * std::sqrt(365.0);
}

BacktestSummary run_backtest(const std::vector<PricePoint>& history, const std::string& model_name,
                              const ModelFactory& factory, int window_days, int trailing_vol_days) {
    if (static_cast<int>(history.size()) < trailing_vol_days + window_days + 1) {
        throw std::runtime_error("run_backtest: not enough history for even one window");
    }

    const double r = 0.0;  // see docs/numerics.md Phase 15: matches how BTC options are
                            // actually quoted (Deribit's own snapshot carries interest_rate: 0.0
                            // on every contract, since collateral is posted in BTC itself).

    std::vector<BacktestWindowResult> windows;

    int window_start_idx = trailing_vol_days;
    const int n = static_cast<int>(history.size());

    while (window_start_idx + window_days < n) {
        double vol = trailing_realized_vol(history, window_start_idx - trailing_vol_days, trailing_vol_days);

        double S0 = history[window_start_idx].price;
        double strike = S0;  // at-the-money at window open
        double T = window_days / 365.0;

        PriceDeltaFn pricer = factory(strike, vol);

        auto [price0, delta0] = pricer(S0, T);

        double cash = -delta0 * S0;
        double shares = delta0;
        const double dt = 1.0 / 365.0;

        for (int step = 1; step < window_days; ++step) {
            int idx = window_start_idx + step;
            double S = history[idx].price;
            double tau = (window_days - step) / 365.0;

            auto [price_t, delta_t] = pricer(S, tau);
            (void)price_t;

            double trade = delta_t - shares;
            cash -= trade * S;
            shares = delta_t;
            cash *= std::exp(r * dt);
        }

        int final_idx = window_start_idx + window_days;
        double S_final = history[final_idx].price;
        cash += shares * S_final;
        cash *= 1.0;  // no growth on the final (already-liquidated) step

        double payoff = std::max(S_final - strike, 0.0);
        double premium_grown = price0 * std::exp(r * T);
        double hedge_pnl = (premium_grown + cash) - payoff;
        double hedge_pnl_bps = hedge_pnl / S0 * 10000.0;

        windows.push_back(BacktestWindowResult{history[window_start_idx].date, S0, strike, price0, S_final,
                                                 payoff, hedge_pnl, hedge_pnl_bps});

        window_start_idx += window_days;
    }

    double sum = 0.0;
    for (const auto& w : windows) sum += w.hedge_pnl_bps;
    double mean_bps = windows.empty() ? 0.0 : sum / windows.size();

    double sq_diff_sum = 0.0;
    for (const auto& w : windows) {
        double d = w.hedge_pnl_bps - mean_bps;
        sq_diff_sum += d * d;
    }
    double stdev_bps = windows.size() > 1 ? std::sqrt(sq_diff_sum / (windows.size() - 1)) : 0.0;

    return BacktestSummary{model_name, static_cast<int>(windows.size()), mean_bps, stdev_bps, windows};
}

}  // namespace deriv
