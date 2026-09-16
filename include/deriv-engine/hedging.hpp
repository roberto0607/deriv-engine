#pragma once

#include "market_data.hpp"
#include "instrument.hpp"
#include <vector>

namespace deriv {

struct HedgeSimResult {
    double final_pnl;
    double theoretical_price;
    std::vector<double> price_path;  // for inspection/plotting
};

// Simulates ONE price path, rebalancing the delta hedge at each of
// `num_rebalances` evenly-spaced points between today and expiry.
HedgeSimResult simulate_delta_hedge(const EuropeanOption& option, const MarketData& market,
                                     int num_rebalances, unsigned int seed = 0);

}  // namespace deriv