#pragma once

#include "market_data.hpp"
#include "instrument.hpp"

namespace deriv {

struct MonteCarloResult {
    double price;
    double standard_error;
};

MonteCarloResult monte_carlo_price(const EuropeanOption& option, const MarketData& market,
                                    int num_paths, bool antithetic = true);

MonteCarloResult monte_carlo_price_control_variate(const EuropeanOption& option, const MarketData& market,
                                                     int num_paths);                                    

}  // namespace deriv