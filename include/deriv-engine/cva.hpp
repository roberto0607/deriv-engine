#pragma once

#include "portfolio.hpp"
#include "bates.hpp"

namespace deriv {

struct CVAResult {
    // Positive = a real cost: how much less this portfolio is worth to
    // you once your counterparty's default risk is priced in.
    double cva;
    double expected_exposure_peak;  // max over time of E[max(portfolio_value, 0)]
    double expected_exposure_avg;   // time-average of the same
    int num_time_steps;
    int num_paths;
};

// Computes a basic, UNILATERAL credit valuation adjustment for
// `portfolio` against a single counterparty:
//
//   CVA = (1 - recovery_rate) * sum_i [ EE(t_i) * discount(t_i) * dPD_i ]
//
// where EE(t_i) = E[max(portfolio_value(t_i), 0)] -- expected POSITIVE
// exposure only, since if the counterparty owes you nothing at time
// t_i, their defaulting costs you nothing; a negative portfolio value
// (you owe THEM) is floored to zero here exactly the way historical_var
// / monte_carlo_var floor nothing but this function specifically must,
// because CVA is about what YOU stand to lose, not the book's raw P&L
// -- discount(t_i) = exp(-market.risk_free_rate * t_i), and dPD_i =
// exp(-hazard_rate*t_{i-1}) - exp(-hazard_rate*t_i) is the marginal
// probability of default landing in that interval under a flat
// (constant) hazard rate.
//
// `num_paths` independent forward paths of the underlying are simulated
// under the Bates SDE (full-truncation Euler for the variance process,
// exactly as in heston_mc.cpp/bates_mc.cpp/var.cpp's own stepping loop
// -- pass jump_intensity=0 in `params` to simulate under pure Heston
// instead) out to `horizon_years`, sampled at `num_time_steps`
// equally-spaced points; the portfolio is repriced -- and its options'
// time-to-expiry decayed via decay_time(), same as everywhere else in
// this project -- at each point on each path to build the expected
// exposure profile EE(t_i).
//
// This is a FOURTH independent implementation of the Bates SDE stepping
// loop (heston_mc.cpp, bates_mc.cpp, and var.cpp's
// simulate_bates_terminal_spot() each already have their own copy) --
// consistent with this project's existing precedent of not sharing
// SDE-stepping code between files (see var.cpp's header comment for the
// same choice), not an oversight.
//
// Deliberately simplified relative to a production CVA calculation --
// see docs/phase-log.md for the honestly-documented list of what's out
// of scope:
//   - Unilateral only: no DVA (your own default risk isn't priced in).
//   - A single flat hazard rate, not a real credit curve bootstrapped
//     from CDS quotes across multiple tenors.
//   - No wrong-way risk: the counterparty's default intensity is
//     modeled as independent of the exposure simulation, when in
//     reality a counterparty's own distress often correlates with
//     exactly the market moves that increase your exposure to them.
//   - One uncollateralized netting set: no CSA/collateral modeling, no
//     netting across multiple counterparties.
//   - No FVA/MVA/KVA -- the other members of the XVA family a real
//     desk computes alongside CVA.
CVAResult compute_cva(const Portfolio& portfolio, const MarketData& market, const BatesParams& params,
                       const PortfolioPricer& pricer, double hazard_rate, double recovery_rate,
                       double horizon_years, int num_time_steps, int num_paths, unsigned int seed = 0);

// Converts a CDS spread (annualized, in decimal -- e.g. 0.02 for 200bps)
// into an approximate flat hazard rate via the standard "credit
// triangle" approximation: hazard_rate ~= spread / (1 - recovery_rate).
// This is itself the simplified version -- a real desk bootstraps a
// full term-structure hazard curve from CDS quotes at several tenors
// rather than reading one flat rate off a single spread, but the
// triangle relationship is the right order-of-magnitude starting point
// and is standard enough to be the textbook first approximation.
double hazard_rate_from_cds_spread(double cds_spread, double recovery_rate);

}  // namespace deriv
