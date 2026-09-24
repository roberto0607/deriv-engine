#pragma once

#include "market_data.hpp"
#include "instrument.hpp"
#include "bates.hpp"
#include "monte_carlo.hpp"

namespace deriv {

// Independent cross-check for bates_price(), mirroring heston_mc.hpp's
// role for heston_price(): simulates the actual Bates SDE (Heston's
// variance process, full-truncation Euler, exactly as in heston_mc.cpp)
// with a compound Poisson jump layered onto the log-price increment each
// step, instead of evaluating the closed-form characteristic function and
// COS series in bates.cpp. This is a genuinely different numerical
// method from the CF/COS path -- it shares no code with
// bates_log_return_cf() or bates_put_price() -- so agreement between the
// two is real evidence the CF derivation (in particular the jump
// compensator k and the compound-Poisson CF identity) matches the model
// it claims to price.
//
// Per step: draw a Poisson(lambda*dt) jump count, sum that many i.i.d.
// N(jump_mean, jump_vol^2) log-jump sizes, and add the total directly to
// the log-price increment alongside the usual diffusion term -- the
// standard Euler discretization of a jump-diffusion SDE. The variance
// process is stepped exactly as in heston_mc.cpp (full truncation); jumps
// don't touch variance in the Bates model, only the log-price.
MonteCarloResult bates_monte_carlo_price(const EuropeanOption& option, const MarketData& market,
                                          const BatesParams& params, int num_paths, int num_steps,
                                          bool antithetic = true, unsigned int seed = 0);

}  // namespace deriv
