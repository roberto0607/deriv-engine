#pragma once

#include "portfolio.hpp"

namespace deriv {

// Portfolio-level gamma/vega/theta via finite difference -- the same
// bump-and-reprice pattern portfolio_delta() already uses in
// portfolio.hpp, extended to the other three Greeks P&L attribution
// needs. Central difference for gamma (second derivative in spot) and
// vega (bump market.volatility); a forward difference for theta (decay
// the portfolio by one calendar day and reprice, since "one day
// earlier" isn't a meaningful bump). Works unchanged for Black-Scholes,
// Heston, or Bates -- exactly like portfolio_delta() and every VaR
// function in var.hpp -- since it only ever calls `pricer` and never
// assumes a closed form.
double portfolio_gamma(const Portfolio& portfolio, const MarketData& market, const PortfolioPricer& pricer);
double portfolio_vega(const Portfolio& portfolio, const MarketData& market, const PortfolioPricer& pricer);
// Sign convention matches greeks.hpp/monte_carlo.cpp: theta is the
// sensitivity to CALENDAR TIME passing (time_to_expiry decreasing), so
// theta = -d(price)/d(time_to_expiry) = +d(price)/d(calendar time). A
// long option's theta is typically negative (it loses value as time
// passes, all else equal) -- same sign a desk would read off a risk
// report.
double portfolio_theta(const Portfolio& portfolio, const MarketData& market, const PortfolioPricer& pricer);

struct PnLAttributionResult {
    double total_pnl;   // actual repriced P&L: value(market_after) - value(market_before)
    double delta_pnl;    // delta * dS
    double gamma_pnl;    // 0.5 * gamma * dS^2
    double vega_pnl;     // vega * dVol (vol POINTS, e.g. 0.60 -> 0.65 is dVol = 0.05, not "5%")
    double theta_pnl;    // theta * dt_years
    // What the four bucket above couldn't account for: cross terms
    // (e.g. spot moving AND vol changing together -- "vanna" -- isn't
    // captured by summing separate delta/vega bumps), higher-order
    // terms the Taylor expansion drops, and any genuine model
    // mis-hedge. A real desk's P&L-explain process treats a large
    // unexplained_pnl as a red flag that the hedge or the model is
    // missing something -- not something to paper over.
    double unexplained_pnl;

    // The portfolio-level Greeks actually used above, computed once at
    // market_before -- exposed for reporting/debugging, the same way a
    // desk's P&L-explain report shows both the Greeks used and the
    // resulting attribution.
    double delta;
    double gamma;
    double vega;
    double theta;
};

// Taylor-decomposes a portfolio's REALIZED P&L between two market
// snapshots into delta/gamma/vega/theta contributions plus an
// unexplained residual -- the standard "P&L explain" process a sell-side
// trading desk runs every single day to check that a position's hedge
// P&L is actually being driven by what its Greeks say it should be
// driven by, not something the model is missing.
//
// Greeks are computed once, at market_before (bump-and-reprice via
// `pricer`, using portfolio_delta()/portfolio_gamma()/portfolio_vega()/
// portfolio_theta() above) -- this is exactly what a desk does in
// practice: yesterday's Greeks are used to predict today's P&L, not
// some Greeks re-derived after the fact from today's market. Because
// this only ever calls `pricer`, it works unchanged for Black-Scholes,
// Heston, or Bates, the same way portfolio_delta() and every function in
// var.hpp already do.
//
// time_elapsed_days must be the actual number of calendar days between
// market_before and market_after (matching how var.hpp's stress tests
// and historical VaR look up real elapsed time from price history)
// -- it decays every option position's time-to-expiry via
// decay_time() before repricing under market_after, so total_pnl
// reflects a real day (or window) passing, not just an instantaneous
// market shock.
PnLAttributionResult attribute_pnl(const Portfolio& portfolio, const MarketData& market_before,
                                    const MarketData& market_after, const PortfolioPricer& pricer,
                                    double time_elapsed_days);

}  // namespace deriv
