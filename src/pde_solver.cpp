#include "deriv-engine/pde_solver.hpp"
#include <cmath>
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

}  // namespace

std::vector<double> solve_tridiagonal(const std::vector<double>& lower,
                                       const std::vector<double>& diag,
                                       const std::vector<double>& upper,
                                       const std::vector<double>& d) {
    std::size_t n = diag.size();
    std::vector<double> c_prime(n);
    std::vector<double> d_prime(n);
    std::vector<double> x(n);

    c_prime[0] = upper[0] / diag[0];
    d_prime[0] = d[0] / diag[0];

    for (std::size_t i = 1; i < n; ++i) {
        double denom = diag[i] - lower[i] * c_prime[i - 1];
        c_prime[i] = upper[i] / denom;
        d_prime[i] = (d[i] - lower[i] * d_prime[i - 1]) / denom;
    }

    x[n - 1] = d_prime[n - 1];
    for (std::size_t i = n - 1; i-- > 0;) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    return x;
}

double crank_nicolson_price(const EuropeanOption& option, const MarketData& market,
                             int num_price_steps, int num_time_steps) {
    const double K = option.strike;
    const double T = option.time_to_expiry;
    const double r = market.risk_free_rate;
    const double q = market.dividend_yield;
    const double sigma = market.volatility;

    const double S_max = 3.0 * K;
    const double dS = S_max / num_price_steps;
    const double dt = T / num_time_steps;

    std::vector<double> V(num_price_steps + 1);
    for (int i = 0; i <= num_price_steps; ++i) {
        double S = i * dS;
        V[i] = payoff(S, K, option.type);
    }

    std::vector<double> alpha(num_price_steps + 1), beta(num_price_steps + 1), gamma(num_price_steps + 1);
    for (int i = 1; i < num_price_steps; ++i) {
        double Si = i * dS;
        double sigma2S2 = sigma * sigma * Si * Si;
        alpha[i] = 0.25 * dt * (sigma2S2 / (dS * dS) - (r - q) * Si / dS);
        beta[i]  = -0.5 * dt * (sigma2S2 / (dS * dS) + r);
        gamma[i] = 0.25 * dt * (sigma2S2 / (dS * dS) + (r - q) * Si / dS);
    }

    for (int step = 0; step < num_time_steps; ++step) {
        std::vector<double> lower(num_price_steps + 1), diag(num_price_steps + 1);
        std::vector<double> upper(num_price_steps + 1), rhs(num_price_steps + 1);

        diag[0] = 1.0;
        upper[0] = 0.0;
        rhs[0] = (option.type == OptionType::Put) ? K * std::exp(-r * (step + 1) * dt) : 0.0;

        diag[num_price_steps] = 1.0;
        lower[num_price_steps] = 0.0;
        rhs[num_price_steps] = (option.type == OptionType::Call)
            ? S_max - K * std::exp(-r * (step + 1) * dt)
            : 0.0;

        for (int i = 1; i < num_price_steps; ++i) {
            lower[i] = -alpha[i];
            diag[i] = 1.0 - beta[i];
            upper[i] = -gamma[i];
            rhs[i] = alpha[i] * V[i - 1] + (1.0 + beta[i]) * V[i] + gamma[i] * V[i + 1];
        }

        V = solve_tridiagonal(lower, diag, upper, rhs);
    }

    double S0 = market.spot;
    int idx = static_cast<int>(S0 / dS);
    idx = std::min(idx, num_price_steps - 1);
    double weight = (S0 - idx * dS) / dS;
    return V[idx] * (1.0 - weight) + V[idx + 1] * weight;
}

double explicit_fd_price(const EuropeanOption& option, const MarketData& market,
                          int num_price_steps, int num_time_steps) {
    const double K = option.strike;
    const double T = option.time_to_expiry;
    const double r = market.risk_free_rate;
    const double q = market.dividend_yield;
    const double sigma = market.volatility;

    const double S_max = 3.0 * K;
    const double dS = S_max / num_price_steps;
    const double dt = T / num_time_steps;

    std::vector<double> V(num_price_steps + 1);
    for (int i = 0; i <= num_price_steps; ++i) {
        V[i] = payoff(i * dS, K, option.type);
    }

    for (int step = 0; step < num_time_steps; ++step) {
        std::vector<double> V_new(num_price_steps + 1);
        V_new[0] = (option.type == OptionType::Put) ? K * std::exp(-r * (step + 1) * dt) : 0.0;
        V_new[num_price_steps] = (option.type == OptionType::Call)
            ? S_max - K * std::exp(-r * (step + 1) * dt)
            : 0.0;

        for (int i = 1; i < num_price_steps; ++i) {
            double Si = i * dS;
            double sigma2S2 = sigma * sigma * Si * Si;
            double a = 0.5 * dt * (sigma2S2 / (dS * dS) - (r - q) * Si / dS);
            double b = 1.0 - dt * (sigma2S2 / (dS * dS) + r);
            double c = 0.5 * dt * (sigma2S2 / (dS * dS) + (r - q) * Si / dS);

            V_new[i] = a * V[i - 1] + b * V[i] + c * V[i + 1];
        }

        V = V_new;
    }

    double S0 = market.spot;
    int idx = static_cast<int>(S0 / dS);
    idx = std::min(idx, num_price_steps - 1);
    double weight = (S0 - idx * dS) / dS;
    return V[idx] * (1.0 - weight) + V[idx + 1] * weight;
}

}  // namespace deriv