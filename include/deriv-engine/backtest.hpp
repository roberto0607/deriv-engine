#pragma once

#include "price_history.hpp"
#include <functional>
#include <string>
#include <vector>

namespace deriv {

// Historical delta-hedging backtest: rolls a fresh at-the-money European
// call every `window_days` calendar days over a real historical price
// series (no simulated paths anywhere in this file), delta-hedges it
// daily using the real historical price path, and records the realized
// hedge P&L at expiry. This is the same accounting hedging.hpp's
// simulate_delta_hedge() uses (collect the premium upfront, track a
// self-financing hedge account, compare total position against the
// option's actual payoff at expiry) -- the only thing that changes here
// is where the price path comes from: real BTC history instead of a
// simulated GBM draw.
//
// Why this needed a new, model-agnostic loop instead of three copies of
// simulate_delta_hedge(): Black-Scholes has a closed-form analytical
// delta, but Heston and Bates don't -- this project's COS-method pricers
// were never extended with an analytical Greeks derivation (a real
// project of its own; see docs/numerics.md Phase 15 for why finite
// differences were used instead). A single generic loop, parameterized
// by a (price, delta) function per model, avoids writing the day-by-day
// rebalancing/cash-accounting logic three times and risking it drifting
// out of sync between models -- exactly the kind of duplication this
// project has avoided elsewhere (see cos_pricing_internal.hpp's shared
// Heston CF for the same reasoning applied to Bates).

// A model's pricer at a single point in time: given the spot price and
// the remaining time to expiry (in years), returns {price, delta}. Each
// BacktestModel below captures whatever it needs (volatility state,
// Heston/Bates structural parameters) as a closure, produced once per
// rolling window by that model's "make a pricer for this window"
// function -- see BacktestWindow.
using PriceDeltaFn = std::function<std::pair<double, double>(double spot, double tau)>;

// Produces a PriceDeltaFn for one rolling window, given that window's
// strike (always set at-the-money: the window's starting spot) and its
// trailing-realized-volatility estimate (annualized, computed from the
// `trailing_vol_days` calendar days immediately before the window
// starts). Each model uses this volatility estimate differently:
// Black-Scholes takes it directly as sigma; Heston/Bates take its square
// as v0 (this window's instantaneous variance) while holding their other
// structural parameters (kappa/theta/xi/rho, and for Bates the jump
// parameters) fixed across the whole backtest -- see backtest.cpp's
// model constructors and docs/numerics.md Phase 15 for why those stay
// fixed rather than being refit per window (there's no historical
// implied-vol surface to refit them against).
using ModelFactory = std::function<PriceDeltaFn(double strike, double window_vol_annualized)>;

struct BacktestWindowResult {
    std::string date_start;
    double spot_start;
    double strike;
    double theoretical_price;
    double spot_end;
    double payoff;
    double hedge_pnl;
    // hedge_pnl expressed in basis points of the window's starting spot
    // (hedge_pnl / spot_start * 10000) -- the metric actually compared
    // and averaged across windows, since BTC's own price spans roughly
    // four orders of magnitude (cents to tens of thousands of dollars)
    // across this dataset's history, making raw dollar P&L meaningless
    // to average across windows.
    double hedge_pnl_bps;
};

struct BacktestSummary {
    std::string model_name;
    int num_windows;
    double mean_pnl_bps;
    double stdev_pnl_bps;
    std::vector<BacktestWindowResult> windows;
};

// Runs the rolling-window backtest described above for one model (via
// its ModelFactory) over the given price history. `window_days` is both
// the option's tenor and the rebalancing interval (daily, at the
// history's own granularity -- see price_history.hpp). `trailing_vol_days`
// is the lookback window for the realized-volatility estimate each new
// window's option is written with. Windows are non-overlapping and
// consume the series strictly in order; the first usable window starts
// once `trailing_vol_days` days of history exist to estimate volatility
// from, and any trailing partial window (fewer than `window_days` days
// remaining) is dropped rather than padded.
BacktestSummary run_backtest(const std::vector<PricePoint>& history, const std::string& model_name,
                              const ModelFactory& factory, int window_days = 30, int trailing_vol_days = 30);

// Annualized realized volatility from the daily log returns of
// history[start_idx .. start_idx + num_days] inclusive of both
// endpoints (so num_days log-returns, from num_days+1 price points),
// using the sample (n-1) standard deviation and a 365-day annualization
// factor (BTC trades every calendar day, unlike traditional equities'
// ~252-day convention). Exposed separately from run_backtest() so it can
// be tested directly against a hand-computed value on a small series.
double trailing_realized_vol(const std::vector<PricePoint>& history, int start_idx, int num_days);

}  // namespace deriv
