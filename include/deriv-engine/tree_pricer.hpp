#pragma once

#include "market_data.hpp"
#include "instrument.hpp"

namespace deriv {

double binomial_tree_price(const EuropeanOption& option, const MarketData& market, int steps);
double binomial_tree_price(const AmericanOption& option, const MarketData& market, int steps);

}  // namespace deriv