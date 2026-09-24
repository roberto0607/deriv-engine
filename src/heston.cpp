#include "deriv-engine/heston.hpp"
#include "deriv-engine/cos_pricing_internal.hpp"
#include <algorithm>
#include <cmath>
#include <complex>

namespace deriv {

namespace {

using cplx = std::complex<double>;
using detail::chi;
using detail::psi;
using detail::heston_log_return_cf;

// COS-method European PUT price under Heston. Calls are recovered via
// put-call parity in heston_price() below, rather than pricing calls
// directly -- see the note on Vk below for why that direction is
// deliberate, not arbitrary.
double heston_put_price(double spot, double strike, double tau, double r, double q,
                         const HestonParams& params) {
    if (tau <= 0.0) {
        return std::max(strike - spot, 0.0);
    }

    // y = ln(S_T / K); the state variable the COS series is expanded
    // in. x = ln(S0 / K) is where "today" sits in that coordinate.
    double x = std::log(spot / strike);

    // Truncation range [a, b] for the COS series. Deliberately generous
    // rather than tightly optimized via exact closed-form Heston
    // cumulants (which are easy to mis-transcribe from memory) -- the
    // COS series converges as long as [a, b] captures essentially all
    // of the terminal log-return distribution's mass, and erring wide
    // only costs a few more (cheap) summation terms, never correctness.
    double avg_var = 0.5 * (params.v0 + params.theta);
    double var_of_var_pad = params.xi * std::sqrt(std::max(avg_var, 1e-10)) * std::sqrt(tau);
    double total_var_estimate = (avg_var + var_of_var_pad) * tau * 2.0;
    double drift_estimate = (r - q - 0.5 * avg_var) * tau;

    const double L = 14.0;
    double width = L * std::sqrt(std::max(total_var_estimate, 1e-10));
    double a = x + drift_estimate - width;
    double b = x + drift_estimate + width;

    // Safety net: the put payoff coefficients integrate over [a, 0],
    // which requires 0 to actually sit inside [a, b]. Guarantee that
    // regardless of how centering/padding worked out above.
    a = std::min(a, -1e-6);
    b = std::max(b, 1e-6);

    const int N = 256;
    double sum = 0.0;

    for (int k = 0; k < N; ++k) {
        double u = k * M_PI / (b - a);
        cplx phi = heston_log_return_cf(u, tau, r, q, params);
        // phi is the CF of the log-return ln(S_T/S0); shifting by x
        // (see derivation notes in heston.hpp's header comment / PR
        // description) turns it into the CF of y = ln(S_T/K) evaluated
        // against this series' basis, without needing a separate
        // "shifted" characteristic function.
        cplx shifted = phi * std::exp(cplx(0.0, 1.0) * u * (x - a));

        double k_pi_over_ba = u;  // = k*pi/(b-a), same quantity

        // Put payoff coefficients integrate over [a, 0] (the region where
        // a put pays off), not [0, b] like a call. This is the fix for a
        // real numerical bug found via property testing: the call-payoff
        // chi/psi evaluate exp(upper) with upper = b, and b grows with
        // total variance (sqrt(vol^2 * tau)) -- for BTC-realistic
        // vol/tenor combinations (e.g. 80% vol, 90d; or 60% vol, 3yr) b
        // routinely lands in the 8-20 range, making exp(b) astronomically
        // large. Each Vk is then O(exp(b)) and the N-term sum has to
        // cancel that down to an O(1)-to-O(spot) answer -- catastrophic
        // cancellation that produced errors from -36% to +∞ (observed:
        // silently wrong prices, not a crash) depending on L and vol/tau,
        // confirmed independent of the summation length N (ruling out
        // under-convergence as the cause). The put's payoff coefficients
        // integrate over [a, 0] instead: exp(0) = 1 and exp(a) is tiny
        // (a is very negative), so Vk stays O(strike) regardless of how
        // wide the truncation range is. Pricing puts directly and calls
        // via put-call parity (below) sidesteps the unstable direction
        // entirely rather than papering over it with a narrower L, which
        // only trades this failure mode for truncation error instead of
        // fixing it.
        double Vk = (2.0 / (b - a)) * strike * (-chi(k_pi_over_ba, a, 0.0, a) + psi(k, k_pi_over_ba, a, 0.0, a));

        double weight = (k == 0) ? 0.5 : 1.0;
        sum += weight * shifted.real() * Vk;
    }

    double price = std::exp(-r * tau) * sum;
    return std::max(price, 0.0);
}

}  // namespace

double heston_price(const EuropeanOption& option, const MarketData& market, const HestonParams& params) {
    double put_price = heston_put_price(market.spot, option.strike, option.time_to_expiry,
                                         market.risk_free_rate, market.dividend_yield, params);

    if (option.type == OptionType::Put) {
        return put_price;
    }

    // Put-call parity: C = P + S*e^{-q*tau} - K*e^{-r*tau}.
    double disc_spot = market.spot * std::exp(-market.dividend_yield * option.time_to_expiry);
    double disc_strike = option.strike * std::exp(-market.risk_free_rate * option.time_to_expiry);
    return put_price + disc_spot - disc_strike;
}

}  // namespace deriv
