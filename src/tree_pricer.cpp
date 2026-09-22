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

// Shared trinomial tree-walking logic, mirroring run_tree() above.
// Boyle (1986): each node branches to three nodes at the next step
// instead of two. The up/down step size dx and the three branch
// probabilities (pu, pm, pd) are chosen so the discretized process
// matches the true lognormal process's mean and variance over one
// time step dt -- the same moment-matching idea behind the binomial
// tree's u/d/p, just with one more free parameter (the middle branch)
// to match with.
//
// Nodes at tree-depth `step` are indexed by j in [-step, step] (2*step+1
// of them), representing price S*exp(j*dx); stored in a vector at
// position i = j + step so indices stay non-negative.
double run_trinomial_tree(double S, double K, double T, double r, double q, double sigma,
                           OptionType type, int steps, bool early_exercise) {
    const double dt = T / steps;
    const double dx = sigma * std::sqrt(3.0 * dt);
    const double nu = r - q - 0.5 * sigma * sigma;
    const double disc = std::exp(-r * dt);

    const double pu = 0.5 * ((sigma * sigma * dt + nu * nu * dt * dt) / (dx * dx) + nu * dt / dx);
    const double pd = 0.5 * ((sigma * sigma * dt + nu * nu * dt * dt) / (dx * dx) - nu * dt / dx);
    const double pm = 1.0 - pu - pd;

    // Prices at expiry: S * exp(j*dx) for j = -steps..steps.
    std::vector<double> values(2 * steps + 1);
    for (int i = 0; i <= 2 * steps; ++i) {
        int j = i - steps;
        double S_final = S * std::exp(j * dx);
        values[i] = payoff(S_final, K, type);
    }

    // Walk backward from expiry to today. At tree-depth `step`, node j
    // (index j+step in the current-level vector) connects to nodes
    // j-1, j, j+1 at depth step+1 (indices j+step, j+step+1, j+step+2
    // in the *previous* iteration's, one-larger, vector).
    for (int step = steps - 1; step >= 0; --step) {
        std::vector<double> next_values(2 * step + 1);
        for (int i = 0; i <= 2 * step; ++i) {
            int j = i - step;
            double continuation = disc * (pu * values[j + step + 2] + pm * values[j + step + 1] +
                                           pd * values[j + step]);

            if (early_exercise) {
                double S_node = S * std::exp(j * dx);
                double exercise_value = payoff(S_node, K, type);
                next_values[i] = std::max(continuation, exercise_value);
            } else {
                next_values[i] = continuation;
            }
        }
        values = std::move(next_values);
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

double trinomial_tree_price(const EuropeanOption& option, const MarketData& market, int steps) {
    return run_trinomial_tree(market.spot, option.strike, option.time_to_expiry,
                               market.risk_free_rate, market.dividend_yield, market.volatility,
                               option.type, steps, /*early_exercise=*/false);
}

double trinomial_tree_price(const AmericanOption& option, const MarketData& market, int steps) {
    return run_trinomial_tree(market.spot, option.strike, option.time_to_expiry,
                               market.risk_free_rate, market.dividend_yield, market.volatility,
                               option.type, steps, /*early_exercise=*/true);
}

}  // namespace deriv