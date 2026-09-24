#pragma once

#include "portfolio.hpp"
#include "price_history.hpp"
#include "bates.hpp"
#include <string>
#include <vector>

namespace deriv {

struct VaRResult {
    std::string method_name;
    int num_scenarios;
    double mean_pnl;
    // Both expressed as POSITIVE numbers -- "the loss at this confidence
    // level" -- even though the underlying P&L distribution is signed.
    // var_95 is exceeded (a worse loss occurs) on average 5% of the
    // time; var_99, 1% of the time.
    double var_95;
    double var_99;
};

// Historical simulation VaR: computes portfolio P&L under every REAL
// horizon_days-ahead return present in `history` (overlapping windows --
// every possible start index, not the non-overlapping windows Phase
// 15's backtest used, since VaR wants the fullest possible empirical
// sample of past moves, not independent draws), applies that same
// percentage move to today's `market.spot`, decrements every option's
// time-to-expiry by horizon_days (via decay_time()) to reflect that time
// has actually passed, and reports the empirical quantiles of the
// resulting P&L. No distributional assumption anywhere in this function
// -- realized history IS the model.
VaRResult historical_var(const Portfolio& portfolio, const MarketData& market,
                          const std::vector<PricePoint>& history, const PortfolioPricer& pricer,
                          int horizon_days = 1);

// Monte Carlo VaR: simulates `num_paths` independent horizon_days-ahead
// terminal spots under the Bates SDE (pass jump_intensity=0 in `params`
// to simulate under pure Heston instead -- the same exact-collapse
// property bates_price() itself has, see tests/test_main.cpp's "[bates]"
// tests), reprices the portfolio at each simulated spot (again via
// decay_time()), and reports the empirical quantiles of simulated P&L.
// Shares no code with historical_var() above and makes a real modeling
// assumption instead of resampling actual history, so agreement between
// the two is genuine cross-validation -- see docs/numerics.md Phase 16.
VaRResult monte_carlo_var(const Portfolio& portfolio, const MarketData& market, const BatesParams& params,
                           const PortfolioPricer& pricer, int num_paths, int horizon_days = 1,
                           unsigned int seed = 0);

// Delta-normal (parametric) VaR: the classic closed-form approximation
// VaR = z(confidence) * |portfolio_delta| * spot * annualized_vol *
// sqrt(horizon_days/365) -- no repricing at all, just the portfolio's
// linear (delta) sensitivity times an assumed-normal move size. z is
// hardcoded for the two confidence levels this project actually uses
// (1.644854 for 95%, 2.326348 for 99% -- standard one-sided normal
// quantiles) rather than a general inverse-normal-CDF solver, since
// nothing else in this codebase needs one yet. Included specifically to
// demonstrate its known blind spot on a book with real convexity (gamma)
// -- see docs/numerics.md Phase 16 -- not as a third method claimed to
// be as trustworthy as the two repricing-based ones above.
double delta_normal_var(const Portfolio& portfolio, const MarketData& market, const PortfolioPricer& pricer,
                         double annualized_vol, int horizon_days, double confidence);

struct StressScenario {
    std::string label;
    std::string start_date;
    std::string end_date;
};

struct StressResult {
    std::string label;
    std::string start_date;
    std::string end_date;
    double spot_start;
    double spot_end;
    double spot_pct_change;
    double portfolio_value_start;
    double portfolio_value_end;
    double pnl;
};

// Replays the portfolio through a real historical window looked up
// directly in `history`: shocks market.spot by the window's ACTUAL
// percentage move, decrements every option's time-to-expiry by the
// actual number of calendar days between start_date and end_date, and
// reports the resulting P&L. No simulation, no distributional
// assumption -- this is a real thing that happened, being replayed
// against today's book.
StressResult run_stress_test(const Portfolio& portfolio, const MarketData& market,
                              const std::vector<PricePoint>& history, const PortfolioPricer& pricer,
                              const StressScenario& scenario);

// Three real, named crash windows already present in
// data/btc_price_history.csv, chosen to cover three different "shapes"
// of loss (see docs/phase-log.md Phase 16): a sharp externally-triggered
// crash (COVID, Feb-Mar 2020), a sharp single-catalyst crash (the FTX
// collapse, Nov 2022), and a long grinding peak-to-trough bear market
// (Nov 2021 top through the Nov 2022 FTX bottom).
std::vector<StressScenario> default_stress_scenarios();

}  // namespace deriv
