#pragma once

#include "market_data.hpp"
#include "instrument.hpp"

namespace deriv {

struct ImpliedVolResult {
    double vol;
    bool converged;
    int iterations;
};

ImpliedVolResult implied_volatility(const EuropeanOption& option, const MarketData& market_without_vol,
                                     double market_price);

}  // namespace deriv