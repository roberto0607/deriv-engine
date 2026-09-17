#include "deriv-engine/implied_vol.hpp"
#include "deriv-engine/black_scholes.hpp"
#include "deriv-engine/greeks.hpp"
#include <cmath>

namespace deriv {

ImpliedVolResult implied_volatility(const EuropeanOption& option, const MarketData& market_without_vol,
                                     double market_price) {
    const double tolerance = 1e-6;
    const int max_iterations = 50;

    // A reasonable starting guess. 50% is a deliberately generic
    // starting point — BTC's typical vol range is wide, but
    // Newton-Raphson converges fast enough that the exact starting
    // value rarely matters much for well-behaved cases.
    double sigma = 0.5;

    for (int i = 0; i < max_iterations; ++i) {
        MarketData trial_market = market_without_vol;
        trial_market.volatility = sigma;

        double price = black_scholes_price(option, trial_market);
        double diff = price - market_price;

        if (std::abs(diff) < tolerance) {
            return ImpliedVolResult{sigma, true, i};
        }

        Greeks g = black_scholes_greeks(option, trial_market);
        double vega = g.vega;

        // If vega is essentially zero (deep ITM/OTM, or near expiry),
        // the Newton-Raphson step is undefined (dividing by ~0) —
        // bail out rather than risk a wild, meaningless jump.
        if (std::abs(vega) < 1e-8) {
            return ImpliedVolResult{sigma, false, i};
        }

        // The actual Newton-Raphson update: adjust sigma by how far
        // off the price is, scaled by how sensitive price is to sigma
        // (vega). This is the same idea as gradient descent, just for
        // a single variable with an exact local derivative available.
        sigma = sigma - diff / vega;

        // Guard against the iteration wandering into nonsensical
        // territory (negative or absurdly large vol) — clamp back into
        // a sane range rather than letting one bad step derail
        // convergence entirely.
        if (sigma <= 0.0) sigma = 0.01;
        if (sigma > 5.0) sigma = 5.0;
    }

    return ImpliedVolResult{sigma, false, max_iterations};
}

}  // namespace deriv