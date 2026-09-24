#pragma once

#include "market_data.hpp"
#include "instrument.hpp"

namespace deriv {

// Path-dependent (Asian, Barrier) option pricing.
//
// Both instruments need the *whole* simulated price path, not just S_T --
// unlike monte_carlo.cpp's vanilla European pricer, which deliberately
// samples only the terminal value (see docs/numerics.md). The path
// simulation itself reuses the same Euler-exact GBM stepping already used
// for American options in longstaff_schwartz.cpp:
//   S[step] = S[step-1] * exp((r-q-0.5*sigma^2)*dt + sigma*sqrt(dt)*z)
//
// Validation strategy (see docs/numerics.md Phase 13 and
// tests/test_main.cpp for the actual tests):
//   - Geometric-average Asian options have a closed form (Kemna-Vorst).
//     There's no closed form for the arithmetic-average Asian option this
//     module actually prices, so the geometric closed form is used as an
//     absolute-correctness anchor for the *path simulation machinery
//     itself* -- if the simulator can reproduce the geometric closed form,
//     the same simulator's arithmetic-average result can be trusted too.
//   - Down-and-out barrier calls (with H <= min(S,K)) have a closed form
//     via the method of images:
//       C_do(S,K,T) = C_BS(S,K,T) - (H/S)^(2*nu) * C_BS(H^2/S, K, T)
//     where nu = (r-q)/sigma^2 - 1/2. Used as an absolute-correctness
//     anchor the same way.
//   - The other three barrier directions (down-and-in, up-and-out,
//     up-and-in) have no closed form implemented here, but the
//     model-independent in/out parity identity
//       (same-direction knock-out) + (same-direction knock-in) == vanilla
//     lets them be validated against the vanilla Black-Scholes price
//     without needing their own formulas.

struct AsianResult {
    double price;
    double standard_error;
};

// Arithmetic-average Asian option via Monte Carlo path simulation.
// avg(S) is the average of num_steps equally-spaced monitoring prices
// from just after today through expiry (inclusive of the terminal price).
AsianResult asian_option_price_mc(const AsianOption& option, const MarketData& market,
                                   int num_paths, int num_steps, unsigned int seed = 0);

// Geometric-average Asian option, closed form (Kemna & Vorst, 1990).
// Exists purely to validate asian_option_price_mc's path-simulation
// machinery -- pass geometric_average=true to asian_option_price_mc and
// the two should agree to within Monte Carlo error.
double geometric_asian_price_closed_form(const AsianOption& option, const MarketData& market,
                                          int num_steps);

// Same Monte Carlo engine as asian_option_price_mc, but averaging with
// the geometric mean instead of the arithmetic mean -- used only to
// validate against geometric_asian_price_closed_form above.
AsianResult geometric_asian_price_mc(const AsianOption& option, const MarketData& market,
                                      int num_paths, int num_steps, unsigned int seed = 0);

struct BarrierResult {
    double price;
    double standard_error;
};

// Barrier option via Monte Carlo path simulation, monitored at num_steps
// equally-spaced dates from just after today through expiry. Handles all
// four BarrierDirection variants.
BarrierResult barrier_option_price_mc(const BarrierOption& option, const MarketData& market,
                                       int num_paths, int num_steps, unsigned int seed = 0);

// Down-and-out call, closed form via the method of images. Only valid
// for a call with H <= min(S, K) -- that constraint is what makes the
// single-reflection image formula exact (see docs/numerics.md Phase 13).
double down_and_out_call_price_closed_form(const BarrierOption& option, const MarketData& market);

}  // namespace deriv
