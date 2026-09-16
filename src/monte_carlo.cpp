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

MonteCarloResult monte_carlo_price_control_variate(const EuropeanOption& option, const MarketData& market,
                                                     int num_paths) {
    const double S = market.spot;
    const double K = option.strike;
    const double T = option.time_to_expiry;
    const double r = market.risk_free_rate;
    const double q = market.dividend_yield;
    const double sigma = market.volatility;

    const double drift = (r - q - 0.5 * sigma * sigma) * T;
    const double diffusion = sigma * std::sqrt(T);

    std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> normal(0.0, 1.0);

    std::vector<double> payoffs(num_paths);
    std::vector<double> controls(num_paths);

    for (int i = 0; i < num_paths; ++i) {
        double z = normal(rng);
        double S_T = S * std::exp(drift + diffusion * z);
        payoffs[i] = payoff(S_T, K, option.type);
        controls[i] = S_T;
    }

    const double control_expected = S * std::exp((r - q) * T);

    double payoff_mean = 0.0, control_mean = 0.0;
    for (int i = 0; i < num_paths; ++i) {
        payoff_mean += payoffs[i];
        control_mean += controls[i];
    }
    payoff_mean /= num_paths;
    control_mean /= num_paths;

    double cov = 0.0, control_var = 0.0;
    for (int i = 0; i < num_paths; ++i) {
        double dp = payoffs[i] - payoff_mean;
        double dc = controls[i] - control_mean;
        cov += dp * dc;
        control_var += dc * dc;
    }
    double c = cov / control_var;

    std::vector<double> adjusted(num_paths);
    for (int i = 0; i < num_paths; ++i) {
        adjusted[i] = payoffs[i] - c * (controls[i] - control_expected);
    }

    double adjusted_mean = 0.0;
    for (double a : adjusted) adjusted_mean += a;
    adjusted_mean /= num_paths;

    double discount = std::exp(-r * T);
    double price = discount * adjusted_mean;

    double sq_diff_sum = 0.0;
    for (double a : adjusted) {
        double diff = a - adjusted_mean;
        sq_diff_sum += diff * diff;
    }
    double sample_variance = sq_diff_sum / (num_paths - 1);
    double standard_error = discount * std::sqrt(sample_variance / num_paths);

    return MonteCarloResult{price, standard_error};
}

}  // namespace derive