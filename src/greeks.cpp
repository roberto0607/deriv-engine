#include "deriv-engine/greeks.hpp"
#include <cmath>

namespace deriv {

namespace {

double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double norm_pdf(double x) {
    static const double inv_sqrt_2pi = 0.3989422804014327;
    return inv_sqrt_2pi * std::exp(-0.5 * x * x);
}

}  // namespace

Greeks black_scholes_greeks(const EuropeanOption& option, const MarketData& market) {
    const double S = market.spot;
    const double K = option.strike;
    const double T = option.time_to_expiry;
    const double r = market.risk_free_rate;
    const double q = market.dividend_yield;
    const double sigma = market.volatility;

    // At or past expiry, there's no more time-sensitivity left: gamma,
    // vega, theta, and rho all collapse to 0. Delta becomes a hard 0/1
    // (or 0/-1 for a put) depending on whether the option finished
    // in-the-money. The formula below would divide by zero
    // (sigma * sqrt(T)) if T <= 0, so this case is handled separately.
    if (T <= 1e-8) {
        Greeks g{};
        const bool call_itm = (option.type == OptionType::Call) && (S > K);
        const bool put_itm  = (option.type == OptionType::Put)  && (S < K);

        if (call_itm) {
            g.delta = 1.0;
        } else if (put_itm) {
            g.delta = -1.0;
        } else {
            g.delta = 0.0;
        }
        // gamma, vega, theta, rho already zero-initialized by Greeks g{};
        return g;
    }

    const double sqrtT = std::sqrt(T);
    const double d1 = (std::log(S / K) + (r - q + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    const double d2 = d1 - sigma * sqrtT;

    const double disc_q = std::exp(-q * T);
    const double disc_r = std::exp(-r * T);
    const double pdf_d1 = norm_pdf(d1);

    Greeks g{};

    // Gamma and vega are identical for calls and puts.
    g.gamma = disc_q * pdf_d1 / (S * sigma * sqrtT);
    g.vega = S * disc_q * pdf_d1 * sqrtT;

    if (option.type == OptionType::Call) {
        g.delta = disc_q * norm_cdf(d1);
        g.theta = -(S * disc_q * pdf_d1 * sigma) / (2.0 * sqrtT)
                  - r * K * disc_r * norm_cdf(d2)
                  + q * S * disc_q * norm_cdf(d1);
        g.rho = K * T * disc_r * norm_cdf(d2);
    } else {
        g.delta = disc_q * (norm_cdf(d1) - 1.0);
        g.theta = -(S * disc_q * pdf_d1 * sigma) / (2.0 * sqrtT)
                  + r * K * disc_r * norm_cdf(-d2)
                  - q * S * disc_q * norm_cdf(-d1);
        g.rho = -K * T * disc_r * norm_cdf(-d2);
    }

    return g;
}

}  // namespace deriv