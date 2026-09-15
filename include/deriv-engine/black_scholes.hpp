#pragma once

#include "market_data.hpp"
#include "instrument.hpp"

namespace deriv {

double black_scholes_price(const EuropeanOption& option, const MarketData& market);

}  // namespace deriv