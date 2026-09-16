#include "deriv-engine/monte_carlo.hpp"
#include "deriv-engine/adouble.hpp"
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

ADouble ad_payoff(const ADouble& S_T, double K, OptionType type) {
    if (type == OptionType::Call) {
        ADouble diff = S_T - K;
        if (diff.value > 0.0) {
            return diff;
        } else {
            return ADouble(0.0);
        }
    } else {
        ADouble diff = ADouble(K) - S_T;
        if (diff.value > 0.0) {
            return diff;
        } else {
            return ADouble(0.0);
        }
    }
}

}  // namespace

MonteCarloResult monte_carlo_price(const EuropeanOption& option, const MarketData& market,
                                    int num_paths, bool antithetic, unsigned int seed) {
    const double S = market.spot;
    const double K = option.strike;
    const double T = option.time_to_expiry;
    const double r = market.risk_free_rate;
    const double q = market.dividend_yield;
    const double sigma = market.volatility;

    const double drift = (r - q - 0.5 * sigma * sigma) * T;
    const double diffusion = sigma * std::sqrt(T);

    // If seed==0 (the default), use a genuinely random seed each call.
    // If a nonzero seed is explicitly passed, use it — this lets two
    // calls (e.g. a base case and a bumped case) reproduce the EXACT
    // same underlying random draws, which finite-difference bumping
    // requires to isolate the effect of the bumped input alone.
    std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
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

MonteCarloGreeksResult monte_carlo_greeks(const EuropeanOption& option, const MarketData& market,
                                           int num_paths) {
    const double K = option.strike;

    std::mt19937 rng(std::random_device{}());
    std::normal_distribution<double> normal(0.0, 1.0);

    double price_sum = 0.0, delta_sum = 0.0, vega_sum = 0.0, theta_sum = 0.0, rho_sum = 0.0;
    std::vector<double> path_prices(num_paths);

    for (int i = 0; i < num_paths; ++i) {
        get_tape().clear();

        ADouble S(market.spot);
        ADouble sigma(market.volatility);
        ADouble r(market.risk_free_rate);
        ADouble T(option.time_to_expiry);
        double q = market.dividend_yield;  // held plain — not a Greek we report here

        double z = normal(rng);  // held FIXED while differentiating — this is what
                                  // makes it "pathwise" differentiation

        ADouble half = ADouble(0.5) * sigma * sigma;
        ADouble drift = (r - ADouble(q) - half) * T;
        ADouble diffusion = sigma * ad_sqrt(T) * z;
        ADouble exponent = drift + diffusion;
        ADouble S_T = S * ad_exp(exponent);

        ADouble path_payoff = ad_payoff(S_T, K, option.type);

        // Discounting is now PART of the tape, not a separate multiply
        // afterward. This means the backward pass automatically captures
        // both ways r affects the price: through the drift (inside the
        // payoff) and through discounting itself — the chain rule handles
        // combining those two effects for us.
        ADouble discount = ad_exp(-(r * T));
        ADouble discounted_price = discount * path_payoff;

        get_tape().backward(discounted_price.idx);
        auto& adj = get_tape().adjoints;

        path_prices[i] = discounted_price.value;
        price_sum += discounted_price.value;
        delta_sum += adj[S.idx];
        vega_sum += adj[sigma.idx];
        rho_sum += adj[r.idx];
        // theta is conventionally the sensitivity to CALENDAR TIME passing,
        // not to time-to-expiry T. Since T = expiry - now, dT/d(calendar
        // time) = -1, so theta = -d(price)/dT.
        theta_sum += -adj[T.idx];
    }

    double price = price_sum / num_paths;
    double delta = delta_sum / num_paths;
    double vega = vega_sum / num_paths;
    double theta = theta_sum / num_paths;
    double rho = rho_sum / num_paths;

    double sq_diff_sum = 0.0;
    for (double p : path_prices) {
        double diff = p - price;
        sq_diff_sum += diff * diff;
    }
    double sample_variance = sq_diff_sum / (num_paths - 1);
    double standard_error = std::sqrt(sample_variance / num_paths);

    return MonteCarloGreeksResult{price, delta, vega, theta, rho, standard_error};
}

}  // namespace derive