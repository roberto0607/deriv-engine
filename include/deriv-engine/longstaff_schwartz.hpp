#pragma once

#include "market_data.hpp"
#include "instrument.hpp"
#include <vector>
#include <array>

namespace deriv {

struct LSMResult {
    double price;
    double standard_error;
};

LSMResult longstaff_schwartz_price(const AmericanOption& option, const MarketData& market,
                                    int num_paths, int num_steps, unsigned int seed = 0);

// Exposed for testing. Fits Y ~ a0 + a1*X + a2*X^2 by least squares.
std::array<double, 3> fit_quadratic(const std::vector<double>& X, const std::vector<double>& Y);

}  // namespace deriv