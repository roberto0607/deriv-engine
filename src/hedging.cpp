#include "deriv-engine/hedging.hpp"
#include "deriv-engine/black_scholes.hpp"
#include "deriv-engine/greeks.hpp"
#include <cmath>
#include <random>

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

HedgeSimResult simulate_delta_hedge(const EuropeanOption& option, const MarketData& market,
                                     int num_rebalances, unsigned int seed) {
    const double K = option.strike;
    const double T = option.time_to_expiry;
    const double r = market.risk_free_rate;
    const double q = market.dividend_yield;
    const double sigma = market.volatility;
    const double dt = T / num_rebalances;

    std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
    std::normal_distribution<double> normal(0.0, 1.0);

    double theoretical_price = black_scholes_price(option, market);

    // Track ONLY the hedging account's cash flow, starting at zero.
    // We'll compare this against the premium separately at the end,
    // rather than mixing the two together throughout the simulation.
    double cash = 0.0;
    double shares_held = 0.0;
    double S = market.spot;
    std::vector<double> price_path;
    price_path.reserve(num_rebalances + 1);
    price_path.push_back(S);

    for (int step = 0; step < num_rebalances; ++step) {
        double time_remaining = T - step * dt;

        EuropeanOption current_option{K, time_remaining, option.type};
        MarketData current_market{S, r, q, sigma};
        Greeks g = black_scholes_greeks(current_option, current_market);

        double target_shares = g.delta;
        double trade = target_shares - shares_held;

        cash -= trade * S;
        shares_held = target_shares;
        cash *= std::exp(r * dt);

        double z = normal(rng);
        double drift = (r - q - 0.5 * sigma * sigma) * dt;
        double diffusion = sigma * std::sqrt(dt) * z;
        S = S * std::exp(drift + diffusion);
        price_path.push_back(S);
    }

    // Liquidate the hedge at expiry.
    cash += shares_held * S;

    // We received `theoretical_price` upfront (grown at the risk-free
    // rate to expiry, since we could have just banked it), and we owe
    // the option's actual payoff. Compare our TOTAL position (premium +
    // hedging account) against that obligation.
    double premium_grown = theoretical_price * std::exp(r * T);
    double obligation = payoff(S, K, option.type);
    double final_pnl = (premium_grown + cash) - obligation;

    return HedgeSimResult{final_pnl, theoretical_price, price_path};
}

}  // namespace deriv