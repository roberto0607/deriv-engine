#pragma once

#include "market_data.hpp"
#include "instrument.hpp"
#include "heston.hpp"

namespace deriv {

// Bates (1996) model: Heston stochastic volatility with a compound
// Poisson jump process (Merton 1976's lognormal jump-size assumption)
// layered on top of the log-price process. Where Heston alone can only
// move the price through continuous diffusion -- however wild the
// vol-of-vol -- Bates adds sudden discontinuous jumps, which is a
// better description of what BTC actually does (a 20% move in an
// afternoon isn't diffusion, however you tune sigma).
//
//   dS/S = (r - q - lambda*k) dt + sqrt(v) dW1 + dJ
//   dv   = kappa (theta - v) dt + xi sqrt(v) dW2,   corr(dW1, dW2) = rho
//   dJ: compound Poisson, intensity lambda (expected jumps/year); each
//       jump multiplies S by (1+Y), ln(1+Y) ~ N(jump_mean, jump_vol^2)
//   k = E[Y] = exp(jump_mean + 0.5*jump_vol^2) - 1  -- the drift
//       compensator that keeps E[S_T] = S0*exp((r-q)*T) exactly, the
//       same no-arbitrage requirement the plain (r-q) drift satisfies
//       in Heston/Black-Scholes without jumps.
//
// Two independent correctness checks, mirroring how heston.hpp is
// validated (see tests/test_main.cpp):
//   1. jump_intensity = 0 must collapse EXACTLY to heston_price() --
//      the jump characteristic function becomes exp(0) = 1 identically.
//   2. In the deterministic-variance limit (xi -> 0, v0 = theta), Bates
//      must collapse to the closed-form Merton (1976) jump-diffusion
//      price -- a Poisson-weighted sum of Black-Scholes prices, derived
//      independently in the test file rather than reusing any of this
//      model's own code.
// A third check (bates_mc.hpp) simulates the actual Bates SDE + jump
// process and must agree with this closed form -- the same
// belt-and-suspenders pattern heston_mc.hpp uses for Heston.
struct BatesParams {
    HestonParams heston;
    double jump_intensity;  // lambda: expected number of jumps per year
    double jump_mean;       // mu_j: mean of ln(1+Y), the log jump size
    double jump_vol;        // sigma_j: std dev of ln(1+Y)
};

// Prices a European option under Bates via the same COS method
// (Fang & Oosterlee, 2008) heston_price() uses, extending Heston's
// characteristic function with the compound Poisson jump term.
double bates_price(const EuropeanOption& option, const MarketData& market, const BatesParams& params);

}  // namespace deriv
