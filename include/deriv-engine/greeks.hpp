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

}  // namespace deriv