#include "deriv-engine/monte_carlo.hpp"
#include <cmath>
#include <random>
#include <vector>

namespace deriv {

namespace {

double payoff(double S, double K, OptionType type) {
    if (type == OptionType::Call) {
        return std::max(S - K, 0.0);
    } else {
        return std::max(K - S, 0.0);
    }
}

}  // namespace

MonteCarloResult monte_carlo_price(const EuropeanOption& option, const MarketData& market,
                                    int num_paths, bool antithetic) {
    const double S = market.spot;
    const double K = option.strike;
    const double T = option.time_to_expiry;
    const double r = market.risk_free_rate;
    const double q = market.dividend_yield;
    const double sigma = market.volatility;

    // For a European option (no path dependency), we don't need to step
    // through many small time increments — GBM's terminal distribution is
    // known exactly, so we can jump directly from today to expiry in one
    // random draw per path. Full path-stepping becomes necessary later
    // for path-dependent payoffs (American exercise via Longstaff-Schwartz).
    const double drift = (r - q - 0.5 * sigma * sigma) * T;
    const double diffusion = sigma * std::sqrt(T);

    std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> normal(0.0, 1.0);

    std::vector<double> payoffs;
    payoffs.reserve(antithetic ? num_paths * 2 : num_paths);

    for (int i = 0; i < num_paths; ++i) {
        double z = normal(rng);

        double S_T = S * std::exp(drift + diffusion * z);
        payoffs.push_back(payoff(S_T, K, option.type));

        if (antithetic) {
            // Reuse the same random draw, negated. This pairs each path
            // with its mirror image, which are negatively correlated —
            // averaging correlated-negative pairs reduces the variance
            // of the overall average, for the same number of random draws.
            double S_T_anti = S * std::exp(drift - diffusion * z);
            payoffs.push_back(payoff(S_T_anti, K, option.type));
        }
    }

    const double discount = std::exp(-r * T);

    double sum = 0.0;
    for (double p : payoffs) sum += p;
    double mean_payoff = sum / payoffs.size();
    double price = discount * mean_payoff;

    double sq_diff_sum = 0.0;
    for (double p : payoffs) {
        double diff = p - mean_payoff;
        sq_diff_sum += diff * diff;
    }
    double sample_variance = sq_diff_sum / (payoffs.size() - 1);
    double standard_error = discount * std::sqrt(sample_variance / payoffs.size());

    return MonteCarloResult{price, standard_error};
}

}  // namespace deriv