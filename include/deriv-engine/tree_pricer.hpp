#pragma once

#include "market_data.hpp"
#include "instrument.hpp"

namespace deriv {

double binomial_tree_price(const EuropeanOption& option, const MarketData& market, int steps);
double binomial_tree_price(const AmericanOption& option, const MarketData& market, int steps);

// Boyle (1986) trinomial tree: same idea as the binomial tree (walk
// backward from expiry, take the discounted expectation at each node),
// but each node branches three ways -- up, middle (unchanged), down --
// instead of two. Provided alongside the binomial tree as an accuracy/
// convergence-rate comparison, not a replacement: see
// docs/phase-log.md Phase 10 for the measured comparison and the
// (deliberate) decision not to make this the default tree.
double trinomial_tree_price(const EuropeanOption& option, const MarketData& market, int steps);
double trinomial_tree_price(const AmericanOption& option, const MarketData& market, int steps);

}  // namespace deriv