#pragma once

#include "market_data.hpp"
#include "instrument.hpp"
#include "heston.hpp"
#include "monte_carlo.hpp"

namespace deriv {

// Independent cross-check for heston_price(): simulates the actual Heston
// SDEs (see heston.hpp) path-by-path via the full-truncation Euler scheme
// (Lord, Koekkoek & van Dijk, 2010) and prices by discounted average
// payoff, instead of evaluating the closed-form characteristic function
// and COS series in heston.cpp. This is deliberately a second, unrelated
// numerical method: if a sign error or dropped term in heston_log_return_cf
// (heston.cpp) silently produced a self-consistent but wrong price, the
// COS method's own tests (which only check the CF against itself, e.g.
// via put-call parity, or against the deterministic-variance BS limit)
// would not catch it -- an actual simulation of the SDEs would.
//
// Full truncation: the variance process is simulated as ordinary Euler
// but every place v enters a square root or the drift's mean-reversion
// term uses v^+ = max(v, 0) instead of v itself. This keeps the
// simulation well-defined even though Euler-discretized CIR variance can
// (and does) go negative between steps -- unlike the exact/QE schemes,
// full truncation is simple enough to audit by eye, at the cost of
// needing more time steps to control discretization bias. That tradeoff
// is fine here: this function exists to validate heston_price(), not to
// replace it, so a converged answer more slowly is enough.
MonteCarloResult heston_monte_carlo_price(const EuropeanOption& option, const MarketData& market,
                                           const HestonParams& params, int num_paths, int num_steps,
                                           bool antithetic = true, unsigned int seed = 0);

}  // namespace deriv
