#include "deriv-engine/exotic_options.hpp"
#include "deriv-engine/black_scholes.hpp"
#include <cmath>
#include <random>
#include <vector>
#include <algorithm>

namespace deriv {

namespace {

double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double payoff(double avg_or_terminal, double K, OptionType type) {
    if (type == OptionType::Call) {
        return std::max(avg_or_terminal - K, 0.0);
    } else {
        return std::max(K - avg_or_terminal, 0.0);
    }
}

// Simulates num_paths independent GBM paths of num_steps equally-spaced
// steps from today to T (paths[p][0] == S0, paths[p][num_steps] == S_T),
// the same Euler-exact stepping longstaff_schwartz.cpp uses.
std::vector<std::vector<double>> simulate_paths(double S0, double r, double q, double sigma,
                                                 double T, int num_paths, int num_steps,
                                                 std::mt19937& rng) {
    std::normal_distribution<double> normal(0.0, 1.0);
    const double dt = T / num_steps;
    const double drift = (r - q - 0.5 * sigma * sigma) * dt;
    const double diffusion = sigma * std::sqrt(dt);

    std::vector<std::vector<double>> paths(num_paths, std::vector<double>(num_steps + 1));
    for (int p = 0; p < num_paths; ++p) {
        paths[p][0] = S0;
        for (int step = 1; step <= num_steps; ++step) {
            double z = normal(rng);
            paths[p][step] = paths[p][step - 1] * std::exp(drift + diffusion * z);
        }
    }
    return paths;
}

double mean_of(const std::vector<double>& xs) {
    double sum = 0.0;
    for (double x : xs) sum += x;
    return sum / xs.size();
}

double standard_error_of(const std::vector<double>& payoffs, double discount) {
    double mean_payoff = mean_of(payoffs);
    double sq_diff_sum = 0.0;
    for (double p : payoffs) {
        double diff = p - mean_payoff;
        sq_diff_sum += diff * diff;
    }
    double sample_variance = sq_diff_sum / (payoffs.size() - 1);
    return discount * std::sqrt(sample_variance / payoffs.size());
}

// Shared engine for both the arithmetic and geometric average Asian
// pricers -- they differ only in how the monitored prices are averaged.
AsianResult asian_price_mc_impl(const AsianOption& option, const MarketData& market,
                                 int num_paths, int num_steps, unsigned int seed,
                                 bool geometric_average) {
    const double S = market.spot;
    const double K = option.strike;
    const double T = option.time_to_expiry;
    const double r = market.risk_free_rate;
    const double q = market.dividend_yield;
    const double sigma = market.volatility;

    std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
    auto paths = simulate_paths(S, r, q, sigma, T, num_paths, num_steps, rng);

    std::vector<double> payoffs;
    payoffs.reserve(num_paths);

    for (int p = 0; p < num_paths; ++p) {
        double avg;
        if (geometric_average) {
            // Average in log-space to avoid needing a running product
            // that could overflow, then exponentiate once at the end --
            // mathematically identical to the geometric mean.
            double log_sum = 0.0;
            for (int step = 1; step <= num_steps; ++step) {
                log_sum += std::log(paths[p][step]);
            }
            avg = std::exp(log_sum / num_steps);
        } else {
            double sum = 0.0;
            for (int step = 1; step <= num_steps; ++step) {
                sum += paths[p][step];
            }
            avg = sum / num_steps;
        }
        payoffs.push_back(payoff(avg, K, option.type));
    }

    const double discount = std::exp(-r * T);
    double price = discount * mean_of(payoffs);
    double standard_error = standard_error_of(payoffs, discount);

    return AsianResult{price, standard_error};
}

}  // namespace

AsianResult asian_option_price_mc(const AsianOption& option, const MarketData& market,
                                   int num_paths, int num_steps, unsigned int seed) {
    return asian_price_mc_impl(option, market, num_paths, num_steps, seed, /*geometric_average=*/false);
}

AsianResult geometric_asian_price_mc(const AsianOption& option, const MarketData& market,
                                      int num_paths, int num_steps, unsigned int seed) {
    return asian_price_mc_impl(option, market, num_paths, num_steps, seed, /*geometric_average=*/true);
}

// Kemna & Vorst (1990): the geometric average of lognormally-distributed
// prices sampled at equally-spaced times is itself lognormal, so the
// geometric Asian option prices exactly like a vanilla Black-Scholes
// option on an underlying with adjusted drift and volatility.
//
//   sigma_hat^2 = sigma^2 * (n+1)(2n+1) / (6*n^2)
//   mu_hat      = (r - q - 0.5*sigma^2)*(n+1)/(2n) + 0.5*sigma_hat^2
//
// where n = num_steps monitoring dates. The adjusted "forward drift"
// mu_hat plays the role of (r_eff - q_eff) in the standard BS formula
// with volatility sigma_hat; discounting is still at the real rate r.
double geometric_asian_price_closed_form(const AsianOption& option, const MarketData& market,
                                          int num_steps) {
    const double S = market.spot;
    const double K = option.strike;
    const double T = option.time_to_expiry;
    const double r = market.risk_free_rate;
    const double q = market.dividend_yield;
    const double sigma = market.volatility;
    const double n = static_cast<double>(num_steps);

    if (T <= 1e-8) {
        return payoff(S, K, option.type);
    }

    const double sigma_hat_sq = sigma * sigma * (n + 1.0) * (2.0 * n + 1.0) / (6.0 * n * n);
    const double sigma_hat = std::sqrt(sigma_hat_sq);
    const double mu_hat = (r - q - 0.5 * sigma * sigma) * (n + 1.0) / (2.0 * n) + 0.5 * sigma_hat_sq;

    const double d1 = (std::log(S / K) + (mu_hat + 0.5 * sigma_hat_sq) * T) / (sigma_hat * std::sqrt(T));
    const double d2 = d1 - sigma_hat * std::sqrt(T);

    // Forward price of the geometric average under the adjusted dynamics,
    // discounted at the real rate r (not mu_hat -- mu_hat only replaces
    // the drift used to build d1/d2, exactly as in the standard BS
    // derivation where r-q sets both the drift and the discount factor's
    // complement; here they're deliberately decoupled).
    const double F = S * std::exp(mu_hat * T);

    if (option.type == OptionType::Call) {
        return std::exp(-r * T) * (F * norm_cdf(d1) - K * norm_cdf(d2));
    } else {
        return std::exp(-r * T) * (K * norm_cdf(-d2) - F * norm_cdf(-d1));
    }
}

BarrierResult barrier_option_price_mc(const BarrierOption& option, const MarketData& market,
                                       int num_paths, int num_steps, unsigned int seed) {
    const double S = market.spot;
    const double K = option.strike;
    const double H = option.barrier;
    const double T = option.time_to_expiry;
    const double r = market.risk_free_rate;
    const double q = market.dividend_yield;
    const double sigma = market.volatility;

    std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
    auto paths = simulate_paths(S, r, q, sigma, T, num_paths, num_steps, rng);

    const bool is_down = (option.direction == BarrierDirection::DownAndOut ||
                          option.direction == BarrierDirection::DownAndIn);
    const bool is_out = (option.direction == BarrierDirection::DownAndOut ||
                         option.direction == BarrierDirection::UpAndOut);

    std::vector<double> payoffs;
    payoffs.reserve(num_paths);

    for (int p = 0; p < num_paths; ++p) {
        bool touched = false;
        // Monitor every simulated step including the starting spot --
        // a barrier at or past today's spot should knock immediately.
        for (int step = 0; step <= num_steps && !touched; ++step) {
            double s = paths[p][step];
            if (is_down ? (s <= H) : (s >= H)) {
                touched = true;
            }
        }

        double S_T = paths[p][num_steps];
        double vanilla_payoff = payoff(S_T, K, option.type);

        // "Out": alive unless touched -> pays the vanilla payoff only if
        // never touched. "In": dead unless touched -> pays the vanilla
        // payoff only if touched at some point.
        bool pays = is_out ? !touched : touched;
        payoffs.push_back(pays ? vanilla_payoff : 0.0);
    }

    const double discount = std::exp(-r * T);
    double price = discount * mean_of(payoffs);
    double standard_error = standard_error_of(payoffs, discount);

    return BarrierResult{price, standard_error};
}

// Method of images for a down-and-out call, valid only for H <= min(S,K)
// (the reflection principle used here relies on the barrier sitting
// below both spot and strike -- above that, the payoff isn't zero at the
// barrier and a single reflection no longer suffices).
//
//   C_do(S,K,T) = C_BS(S,K,T) - (H/S)^(2*nu) * C_BS(H^2/S, K, T)
//   nu = (r-q)/sigma^2 - 1/2
double down_and_out_call_price_closed_form(const BarrierOption& option, const MarketData& market) {
    const double S = market.spot;
    const double K = option.strike;
    const double H = option.barrier;
    const double r = market.risk_free_rate;
    const double q = market.dividend_yield;
    const double sigma = market.volatility;

    EuropeanOption vanilla{option.strike, option.time_to_expiry, OptionType::Call};
    double C_vanilla = black_scholes_price(vanilla, market);

    // Already knocked out: the barrier is at or above spot, so a
    // down-and-out option monitored continuously would already be dead.
    if (H >= S) {
        return 0.0;
    }

    const double nu = (r - q) / (sigma * sigma) - 0.5;

    MarketData image_market = market;
    image_market.spot = (H * H) / S;
    double C_image = black_scholes_price(vanilla, image_market);

    double factor = std::pow(H / S, 2.0 * nu);

    return C_vanilla - factor * C_image;
}

}  // namespace deriv
