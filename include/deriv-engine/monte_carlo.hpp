#pragma once

#include "market_data.hpp"
#include "instrument.hpp"

namespace deriv {

struct MonteCarloResult {
    double price;
    double standard_error;
};

MonteCarloResult monte_carlo_price(const EuropeanOption& option, const MarketData& market,
                                    int num_paths, bool antithetic = true, unsigned int seed = 0);

MonteCarloResult monte_carlo_price_control_variate(const EuropeanOption& option, const MarketData& market,
                                                     int num_paths);

struct MonteCarloGreeksResult {
    double price;
    double delta;
    double vega;
    double theta;
    double rho;
    double standard_error;
};

MonteCarloGreeksResult monte_carlo_greeks(const EuropeanOption& option, const MarketData& market,
                                           int num_paths);

struct MonteCarloGammaResult {
    double price;
    double delta;
    double gamma;
    double standard_error;
};

// Second-order AAD: gamma via forward-over-reverse automatic
// differentiation (see adouble2.hpp/tape2.hpp), not finite differences.
// Delta is recomputed here too (rather than reusing monte_carlo_greeks'
// delta) as a same-tape consistency check -- see tests/test_main.cpp
// and docs/phase-log.md Phase 11 for why that check matters and what
// it's cross-validated against.
MonteCarloGammaResult monte_carlo_gamma(const EuropeanOption& option, const MarketData& market,
                                         int num_paths, unsigned int seed = 0);

}  // namespace deriv