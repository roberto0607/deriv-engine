#include "deriv-engine/cva.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

namespace deriv {

CVAResult compute_cva(const Portfolio& portfolio, const MarketData& market, const BatesParams& params,
                       const PortfolioPricer& pricer, double hazard_rate, double recovery_rate,
                       double horizon_years, int num_time_steps, int num_paths, unsigned int seed) {
    if (num_time_steps <= 0) throw std::runtime_error("compute_cva: num_time_steps must be positive");
    if (num_paths <= 0) throw std::runtime_error("compute_cva: num_paths must be positive");
    if (horizon_years <= 0.0) throw std::runtime_error("compute_cva: horizon_years must be positive");

    const double dt_years = horizon_years / num_time_steps;
    const double dt_days = dt_years * 365.0;
    const auto& hp = params.heston;

    // The portfolio as it will look at each future sample point --
    // decay_time() is measured from today, so these are precomputed
    // once (not per path) exactly the way var.cpp precomputes a single
    // `decayed` portfolio for its one horizon.
    std::vector<Portfolio> decayed_at_step;
    decayed_at_step.reserve(num_time_steps + 1);
    for (int step = 0; step <= num_time_steps; ++step) {
        decayed_at_step.push_back(decay_time(portfolio, dt_days * step));
    }

    // Expected (positive) exposure at each sample point, index 0 = today
    // (deterministic -- no simulation needed for "right now").
    std::vector<double> expected_exposure(num_time_steps + 1, 0.0);
    expected_exposure[0] = std::max(portfolio_value(portfolio, market, pricer), 0.0);

    std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
    std::normal_distribution<double> normal(0.0, 1.0);

    double k = std::exp(params.jump_mean + 0.5 * params.jump_vol * params.jump_vol) - 1.0;
    double r_compensated = market.risk_free_rate - params.jump_intensity * k;
    double sqrt_1_minus_rho_sq = std::sqrt(std::max(1.0 - hp.rho * hp.rho, 0.0));
    std::poisson_distribution<int> poisson(params.jump_intensity * dt_years);
    const double sqrt_dt = std::sqrt(dt_years);

    std::vector<double> exposure_sum(num_time_steps + 1, 0.0);  // index 0 unused

    for (int path = 0; path < num_paths; ++path) {
        double v = hp.v0;
        double log_S = std::log(market.spot);

        for (int step = 1; step <= num_time_steps; ++step) {
            double v_plus = std::max(v, 0.0);
            double sqrt_v_plus = std::sqrt(v_plus);

            double z1 = normal(rng);
            double z2 = normal(rng);
            double z_var = hp.rho * z1 + sqrt_1_minus_rho_sq * z2;

            double jump_sum = 0.0;
            int num_jumps = poisson(rng);
            for (int j = 0; j < num_jumps; ++j) {
                jump_sum += params.jump_mean + params.jump_vol * normal(rng);
            }

            log_S += (r_compensated - market.dividend_yield - 0.5 * v_plus) * dt_years +
                     sqrt_v_plus * sqrt_dt * z1 + jump_sum;
            v = v + hp.kappa * (hp.theta - v_plus) * dt_years + hp.xi * sqrt_v_plus * sqrt_dt * z_var;

            double S_t = std::exp(log_S);
            MarketData m = market;
            m.spot = S_t;

            double value = portfolio_value(decayed_at_step[step], m, pricer);
            exposure_sum[step] += std::max(value, 0.0);
        }
    }

    for (int step = 1; step <= num_time_steps; ++step) {
        expected_exposure[step] = exposure_sum[step] / num_paths;
    }

    double cva_before_lgd = 0.0;
    double peak = expected_exposure[0];
    double sum_ee = expected_exposure[0];

    for (int step = 1; step <= num_time_steps; ++step) {
        double t_prev = dt_years * (step - 1);
        double t_curr = dt_years * step;

        double survival_prev = std::exp(-hazard_rate * t_prev);
        double survival_curr = std::exp(-hazard_rate * t_curr);
        double default_prob_interval = survival_prev - survival_curr;

        double discount = std::exp(-market.risk_free_rate * t_curr);

        cva_before_lgd += expected_exposure[step] * discount * default_prob_interval;

        peak = std::max(peak, expected_exposure[step]);
        sum_ee += expected_exposure[step];
    }

    double cva = (1.0 - recovery_rate) * cva_before_lgd;
    double avg_ee = sum_ee / (num_time_steps + 1);

    return CVAResult{cva, peak, avg_ee, num_time_steps, num_paths};
}

double hazard_rate_from_cds_spread(double cds_spread, double recovery_rate) {
    if (recovery_rate >= 1.0) throw std::runtime_error("hazard_rate_from_cds_spread: recovery_rate must be < 1");
    return cds_spread / (1.0 - recovery_rate);
}

}  // namespace deriv
