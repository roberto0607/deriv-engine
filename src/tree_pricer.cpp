#include "deriv-engine/tree_pricer.hpp"
#include <cmath>
#include <vector>
#include <algorithm>

namespace deriv {

namespace {

double payoff(double S, double K, OptionType type) {
    if (type == OptionType::Call) {
        return std::max(S - K, 0.0);
    } else {
        return std::max(K - S, 0.0);
    }
}

// Shared tree-walking logic. early_exercise=false reproduces European
// pricing; early_exercise=true checks intrinsic value at every node,
// reproducing American pricing.
double run_tree(double S, double K, double T, double r, double q, double sigma,
                 OptionType type, int steps, bool early_exercise) {
    const double dt = T / steps;
    const double u = std::exp(sigma * std::sqrt(dt));
    const double d = 1.0 / u;
    const double disc = std::exp(-r * dt);
    const double p = (std::exp((r - q) * dt) - d) / (u - d);

    // Prices at the final step (expiry): S * u^j * d^(steps-j) for j = 0..steps
    std::vector<double> values(steps + 1);
    for (int j = 0; j <= steps; ++j) {
        double S_final = S * std::pow(u, j) * std::pow(d, steps - j);
        values[j] = payoff(S_final, K, type);
    }

    // Walk backward from expiry to today, one time-step at a time.
    for (int step = steps - 1; step >= 0; --step) {
        for (int j = 0; j <= step; ++j) {
            double continuation = disc * (p * values[j + 1] + (1.0 - p) * values[j]);

            if (early_exercise) {
                double S_node = S * std::pow(u, j) * std::pow(d, step - j);
                double exercise_value = payoff(S_node, K, type);
                values[j] = std::max(continuation, exercise_value);
            } else {
                values[j] = continuation;
            }
        }
    }

    return values[0];
}

}  // namespace

double binomial_tree_price(const EuropeanOption& option, const MarketData& market, int steps) {
    return run_tree(market.spot, option.strike, option.time_to_expiry,
                     market.risk_free_rate, market.dividend_yield, market.volatility,
                     option.type, steps, /*early_exercise=*/false);
}

double binomial_tree_price(const AmericanOption& option, const MarketData& market, int steps) {
    return run_tree(market.spot, option.strike, option.time_to_expiry,
                     market.risk_free_rate, market.dividend_yield, market.volatility,
                     option.type, steps, /*early_exercise=*/true);
}

}  // namespace deriv