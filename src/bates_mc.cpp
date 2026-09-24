#include "deriv-engine/bates_mc.hpp"
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

// One path of (S, v) from t=0 to t=T over num_steps full-truncation Euler
// steps (identical variance-process stepping to heston_mc.cpp's
// simulate_terminal_spot), with a compound Poisson jump added to the
// log-price increment each step. jump_log_sum[step] is the sum of that
// step's jump-size draws (already ln(1+Y_i) values, pre-summed by the
// caller so this function stays a pure "step the SDE" routine, matching
// heston_mc.cpp's structure) -- zero when no jumps landed in that step.
// The drift uses r_compensated = r - lambda*k, the same compensation
// bates_log_return_cf() applies, so the discretized process is consistent
// with the same risk-neutral martingale condition the closed form prices
// under.
double simulate_terminal_spot(double S0, double v0, double r_compensated, double q, const HestonParams& p,
                               double T, int num_steps, const std::vector<double>& z_spot,
                               const std::vector<double>& z_var_indep,
                               const std::vector<double>& jump_log_sum) {
    const double dt = T / num_steps;
    const double sqrt_dt = std::sqrt(dt);
    const double sqrt_1_minus_rho_sq = std::sqrt(std::max(1.0 - p.rho * p.rho, 0.0));

    double S = S0;
    double v = v0;
    double log_S = std::log(S0);

    for (int step = 0; step < num_steps; ++step) {
        double v_plus = std::max(v, 0.0);
        double sqrt_v_plus = std::sqrt(v_plus);

        double z1 = z_spot[step];
        double z2 = z_var_indep[step];
        double z_var = p.rho * z1 + sqrt_1_minus_rho_sq * z2;

        log_S += (r_compensated - q - 0.5 * v_plus) * dt + sqrt_v_plus * sqrt_dt * z1 + jump_log_sum[step];
        double v_next = v + p.kappa * (p.theta - v_plus) * dt + p.xi * sqrt_v_plus * sqrt_dt * z_var;

        v = v_next;
    }

    (void)S;
    return std::exp(log_S);
}

}  // namespace

MonteCarloResult bates_monte_carlo_price(const EuropeanOption& option, const MarketData& market,
                                          const BatesParams& params, int num_paths, int num_steps,
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

    // Same jump compensator k and compensated drift bates_log_return_cf()
    // uses, so the simulated process and the closed form target the same
    // risk-neutral forward.
    double k = std::exp(params.jump_mean + 0.5 * params.jump_vol * params.jump_vol) - 1.0;
    double r_compensated = r - params.jump_intensity * k;

    std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    std::poisson_distribution<int> poisson(params.jump_intensity * (T / num_steps));

    std::vector<double> z_spot(num_steps), z_var(num_steps), jump_log_sum(num_steps);
    std::vector<double> payoffs;
    payoffs.reserve(antithetic ? static_cast<size_t>(num_paths) * 2 : num_paths);

    for (int i = 0; i < num_paths; ++i) {
        for (int step = 0; step < num_steps; ++step) {
            z_spot[step] = normal(rng);
            z_var[step] = normal(rng);

            int num_jumps = poisson(rng);
            double sum = 0.0;
            for (int j = 0; j < num_jumps; ++j) {
                sum += params.jump_mean + params.jump_vol * normal(rng);
            }
            jump_log_sum[step] = sum;
        }

        double S_T = simulate_terminal_spot(S0, params.heston.v0, r_compensated, q, params.heston, T, num_steps,
                                             z_spot, z_var, jump_log_sum);
        payoffs.push_back(payoff(S_T, K, option.type));

        if (antithetic) {
            // As in heston_mc.cpp: negating both diffusion normal
            // sequences gives a valid antithetic path for the continuous
            // part (the correlation transform is linear). The jump
            // process is left as-is -- it isn't Gaussian-symmetric in
            // the same way, and re-using the same jump draws for both
            // legs of the pair still reduces variance on the diffusion
            // component without needing a jump-specific antithetic
            // scheme.
            std::vector<double> z_spot_anti(num_steps), z_var_anti(num_steps);
            for (int step = 0; step < num_steps; ++step) {
                z_spot_anti[step] = -z_spot[step];
                z_var_anti[step] = -z_var[step];
            }
            double S_T_anti = simulate_terminal_spot(S0, params.heston.v0, r_compensated, q, params.heston, T,
                                                       num_steps, z_spot_anti, z_var_anti, jump_log_sum);
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
