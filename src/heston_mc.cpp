#include "deriv-engine/heston_mc.hpp"
#include <algorithm>
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

// Simulates one path of (S, v) from t=0 to t=T over num_steps full-
// truncation Euler steps, given two independent standard-normal draw
// sequences (one per step, for the spot and variance drivers before
// correlating them). Returns the terminal spot price.
double simulate_terminal_spot(double S0, double v0, double r, double q, const HestonParams& p, double T,
                               int num_steps, const std::vector<double>& z_spot,
                               const std::vector<double>& z_var_indep) {
    const double dt = T / num_steps;
    const double sqrt_dt = std::sqrt(dt);
    const double sqrt_1_minus_rho_sq = std::sqrt(std::max(1.0 - p.rho * p.rho, 0.0));

    double S = S0;
    double v = v0;

    for (int step = 0; step < num_steps; ++step) {
        double v_plus = std::max(v, 0.0);
        double sqrt_v_plus = std::sqrt(v_plus);

        double z1 = z_spot[step];
        double z2 = z_var_indep[step];
        double z_var = p.rho * z1 + sqrt_1_minus_rho_sq * z2;  // correlated with z1

        double S_next = S * std::exp((r - q - 0.5 * v_plus) * dt + sqrt_v_plus * sqrt_dt * z1);
        double v_next = v + p.kappa * (p.theta - v_plus) * dt + p.xi * sqrt_v_plus * sqrt_dt * z_var;

        S = S_next;
        v = v_next;
    }

    return S;
}

}  // namespace

MonteCarloResult heston_monte_carlo_price(const EuropeanOption& option, const MarketData& market,
                                           const HestonParams& params, int num_paths, int num_steps,
                                           bool antithetic, unsigned int seed) {
    const double S0 = market.spot;
    const double K = option.strike;
    const double T = option.time_to_expiry;
    const double r = market.risk_free_rate;
    const double q = market.dividend_yield;

    if (T <= 0.0) {
        double intrinsic = payoff(S0, K, option.type);
        return MonteCarloResult{intrinsic, 0.0};
    }

    std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
    std::normal_distribution<double> normal(0.0, 1.0);

    std::vector<double> z_spot(num_steps), z_var(num_steps);
    std::vector<double> payoffs;
    payoffs.reserve(antithetic ? static_cast<size_t>(num_paths) * 2 : num_paths);

    for (int i = 0; i < num_paths; ++i) {
        for (int step = 0; step < num_steps; ++step) {
            z_spot[step] = normal(rng);
            z_var[step] = normal(rng);
        }

        double S_T = simulate_terminal_spot(S0, params.v0, r, q, params, T, num_steps, z_spot, z_var);
        payoffs.push_back(payoff(S_T, K, option.type));

        if (antithetic) {
            // Negating both underlying normal sequences flips the sign of
            // every z1 and z2 the correlation transform above is built
            // from, which flips zS and zV together correctly (the
            // transform is linear), giving a valid antithetic path rather
            // than an inconsistent one.
            std::vector<double> z_spot_anti(num_steps), z_var_anti(num_steps);
            for (int step = 0; step < num_steps; ++step) {
                z_spot_anti[step] = -z_spot[step];
                z_var_anti[step] = -z_var[step];
            }
            double S_T_anti = simulate_terminal_spot(S0, params.v0, r, q, params, T, num_steps, z_spot_anti,
                                                       z_var_anti);
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
