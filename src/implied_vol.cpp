#include "deriv-engine/implied_vol.hpp"
#include "deriv-engine/black_scholes.hpp"
#include "deriv-engine/greeks.hpp"
#include <algorithm>
#include <cmath>

namespace deriv {

namespace {

// Brenner-Subrahmanyam approximation: a closed-form first guess for
// implied vol, derived from the ATM Black-Scholes price expansion
// (price ~= spot * sigma * sqrt(T / 2pi)). It's only exact at-the-money,
// but as a *starting point* for Newton-Raphson it lands far closer to
// the true root than a flat guess across the whole moneyness range,
// which is what cuts down non-convergence on deep ITM/OTM contracts.
double initial_guess(double option_price, double spot, double time_to_expiry) {
    if (spot <= 0.0 || time_to_expiry <= 0.0) return 0.5;
    double guess = option_price / spot * std::sqrt(2.0 * M_PI / time_to_expiry);
    // Clamp into a sane range -- the approximation can misbehave for
    // options that are far from ATM, where it isn't really valid.
    if (!std::isfinite(guess) || guess < 0.01) guess = 0.3;
    if (guess > 5.0) guess = 5.0;
    return guess;
}

// Bisection: slower than Newton-Raphson but can't diverge or get stuck
// on near-zero vega. Used only as a fallback once Newton-Raphson stalls,
// so its cost is paid only on the hard cases (deep ITM/OTM, near expiry)
// that were previously abandoned outright.
ImpliedVolResult bisection_fallback(const EuropeanOption& option, const MarketData& market_without_vol,
                                     double market_price, int iterations_so_far) {
    const double tolerance = 1e-6;
    const int max_iterations = 100;

    double lo = 1e-4;
    double hi = 5.0;

    auto price_at = [&](double sigma) {
        MarketData trial = market_without_vol;
        trial.volatility = sigma;
        return black_scholes_price(option, trial);
    };

    double price_lo = price_at(lo);
    double price_hi = price_at(hi);

    // If the target price isn't bracketed by [lo, hi] at all, no root
    // exists in range -- e.g. a price quote outside what's achievable
    // for any vol (data noise, a stale quote). Report non-convergence
    // rather than returning a meaningless answer.
    if ((price_lo - market_price) * (price_hi - market_price) > 0.0) {
        return ImpliedVolResult{std::clamp(market_price > price_hi ? hi : lo, lo, hi), false,
                                 iterations_so_far};
    }

    // Black-Scholes price is monotonically increasing in sigma (vega >= 0
    // always), so the direction to narrow in is decided directly by
    // comparing price_mid to the target -- no sign-matching trick needed.
    for (int i = 0; i < max_iterations; ++i) {
        double mid = 0.5 * (lo + hi);
        double price_mid = price_at(mid);
        double diff = price_mid - market_price;

        if (std::abs(diff) < tolerance || (hi - lo) < 1e-9) {
            return ImpliedVolResult{mid, true, iterations_so_far + i};
        }

        if (price_mid < market_price) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return ImpliedVolResult{0.5 * (lo + hi), false, iterations_so_far + max_iterations};
}

}  // namespace

ImpliedVolResult implied_volatility(const EuropeanOption& option, const MarketData& market_without_vol,
                                     double market_price) {
    const double tolerance = 1e-6;
    const int max_iterations = 50;

    double sigma = initial_guess(market_price, market_without_vol.spot, option.time_to_expiry);

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
        // the Newton-Raphson step is undefined (dividing by ~0). Rather
        // than bailing out here, hand off to bisection -- it has no
        // derivative to divide by, so it can still make progress on
        // exactly the cases where Newton-Raphson can't.
        if (std::abs(vega) < 1e-8) {
            return bisection_fallback(option, market_without_vol, market_price, i);
        }

        // The actual Newton-Raphson update: adjust sigma by how far
        // off the price is, scaled by how sensitive price is to sigma
        // (vega). This is the same idea as gradient descent, just for
        // a single variable with an exact local derivative available.
        sigma = sigma - diff / vega;

        // Guard against the iteration wandering into nonsensical
        // territory (negative or absurdly large vol) -- clamp back into
        // a sane range rather than letting one bad step derail
        // convergence entirely.
        if (sigma <= 0.0) sigma = 0.01;
        if (sigma > 5.0) sigma = 5.0;
    }

    // Newton-Raphson used its full iteration budget without converging
    // (e.g. oscillating near a kink) -- try bisection once as a last
    // resort before giving up on this contract entirely.
    return bisection_fallback(option, market_without_vol, market_price, max_iterations);
}

}  // namespace deriv
