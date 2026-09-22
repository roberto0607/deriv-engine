#pragma once

#include "market_data.hpp"
#include "instrument.hpp"

namespace deriv {

struct Greeks {
    double delta;
    double gamma;
    double vega;
    double theta;
    double rho;
};

Greeks black_scholes_greeks(const EuropeanOption& option, const MarketData& market);

struct GammaAD2Result {
    double price;
    double delta;
    double gamma;
};

// Price, delta AND gamma, all read off a single second-order AAD
// (forward-over-reverse, see adouble2.hpp/tape2.hpp) backward pass
// through the Black-Scholes closed-form formula itself -- not the
// Monte Carlo simulation monte_carlo_gamma() also computes gamma for.
// Exists specifically as the working half of the gamma-via-AAD story:
// applying the exact same technique to Monte Carlo's kinked per-path
// payoff gives exactly zero (see monte_carlo_gamma(), and
// docs/phase-log.md Phase 11 for why). This function's gamma should
// match black_scholes_greeks()'s analytical gamma to near machine
// precision -- see tests/test_main.cpp.
GammaAD2Result black_scholes_greeks_via_ad2(const EuropeanOption& option, const MarketData& market);

}  // namespace deriv