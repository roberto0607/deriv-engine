#include "deriv-engine/black_scholes.hpp"
#include <cmath>

namespace deriv {

namespace {

double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

}  // namespace

double black_scholes_price(const EuropeanOption& option, const MarketData& market) {
    const double S = market.spot;
    const double K = option.strike;
    const double T = option.time_to_expiry;
    const double r = market.risk_free_rate;
    const double q = market.dividend_yield;
    const double sigma = market.volatility;

    const double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    const double d2 = d1 - sigma * std::sqrt(T);

    if (option.type == OptionType::Call) {
        return S * std::exp(-q * T) * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
    } else {
        return K * std::exp(-r * T) * norm_cdf(-d2) - S * std::exp(-q * T) * norm_cdf(-d1);
    }
}

}  // namespace deriv