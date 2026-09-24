#include "deriv-engine/bates.hpp"
#include "deriv-engine/cos_pricing_internal.hpp"
#include <algorithm>
#include <cmath>
#include <complex>

namespace deriv {

namespace {

using cplx = std::complex<double>;
using detail::chi;
using detail::psi;

// Compound Poisson jump component of the log-return characteristic
// function: E[exp(iu * sum_{i=1}^{N(tau)} X_i)] where N(tau) ~
// Poisson(lambda*tau) and each X_i ~ N(jump_mean, jump_vol^2) i.i.d.
// (the log of one jump multiplier, ln(1+Y_i)). Standard compound-Poisson
// CF identity: exp(lambda*tau*(phi_X(u) - 1)), phi_X(u) = the CF of a
// single normal jump size.
cplx jump_log_return_cf(double u, double tau, const BatesParams& p) {
    const cplx i(0.0, 1.0);
    cplx phi_single_jump = std::exp(i * u * p.jump_mean - 0.5 * u * u * p.jump_vol * p.jump_vol);
    return std::exp(p.jump_intensity * tau * (phi_single_jump - 1.0));
}

// Characteristic function of the log-return ln(S_T/S0) under Bates:
// Heston's CF (reused verbatim from cos_pricing_internal.hpp, not
// re-derived) times the jump CF above, with the drift rate adjusted by
// -lambda*k so the jump process doesn't change the risk-neutral forward
// price. Passing (r - lambda*k) into heston_log_return_cf achieves that
// compensation directly, since r only ever enters that function through
// the iu*(r-q)*tau drift term -- no separate compensator term is needed
// in jump_log_return_cf above.
cplx bates_log_return_cf(double u, double tau, double r, double q, const BatesParams& p) {
    double k = std::exp(p.jump_mean + 0.5 * p.jump_vol * p.jump_vol) - 1.0;
    double r_compensated = r - p.jump_intensity * k;
    cplx heston_part = detail::heston_log_return_cf(u, tau, r_compensated, q, p.heston);
    cplx jump_part = jump_log_return_cf(u, tau, p);
    return heston_part * jump_part;
}

// COS-method European PUT price under Bates. Calls are recovered via
// put-call parity in bates_price() below, for the same reason
// heston_put_price() prices puts directly (see heston.cpp): the put's
// payoff coefficients stay O(strike) regardless of truncation width,
// avoiding the catastrophic cancellation that pricing deep-ITM calls
// directly (or deriving small calls from large puts) can produce.
double bates_put_price(double spot, double strike, double tau, double r, double q,
                        const BatesParams& params) {
    if (tau <= 0.0) {
        return std::max(strike - spot, 0.0);
    }

    double x = std::log(spot / strike);

    // Truncation range, same padding strategy as heston_put_price()
    // (generous rather than tightly optimized), plus an extra term for
    // the jump process's own contribution to the terminal log-return's
    // variance: Var(compound Poisson sum) = lambda*tau*E[X^2] =
    // lambda*tau*(jump_mean^2 + jump_vol^2), the standard second-moment
    // formula for a compound Poisson process. Without this, a
    // meaningfully-sized jump_vol (e.g. 0.3-0.5, realistic for BTC) can
    // put real probability mass of the terminal distribution outside
    // [a,b], truncating it away and silently underpricing far wings.
    double avg_var = 0.5 * (params.heston.v0 + params.heston.theta);
    double var_of_var_pad = params.heston.xi * std::sqrt(std::max(avg_var, 1e-10)) * std::sqrt(tau);
    double jump_var = params.jump_intensity * tau *
                       (params.jump_mean * params.jump_mean + params.jump_vol * params.jump_vol);
    double total_var_estimate = (avg_var + var_of_var_pad) * tau * 2.0 + jump_var;

    // Truncation-window centering only (the CF itself carries the exact
    // drift, this just keeps [a,b] pointed at the right place): the mean
    // log-return is the usual (r - q - lambda*k - 0.5*avg_var)*tau from
    // the compensated diffusion, plus lambda*tau*jump_mean from the
    // jumps' own expected contribution (E[N(tau)] * E[X] = lambda*tau*jump_mean).
    double k = std::exp(params.jump_mean + 0.5 * params.jump_vol * params.jump_vol) - 1.0;
    double drift_estimate = (r - q - params.jump_intensity * k - 0.5 * avg_var) * tau +
                             params.jump_intensity * tau * params.jump_mean;

    const double L = 14.0;
    double width = L * std::sqrt(std::max(total_var_estimate, 1e-10));
    double a = x + drift_estimate - width;
    double b = x + drift_estimate + width;

    a = std::min(a, -1e-6);
    b = std::max(b, 1e-6);

    const int N = 256;
    double sum = 0.0;

    for (int k_idx = 0; k_idx < N; ++k_idx) {
        double u = k_idx * M_PI / (b - a);
        cplx phi = bates_log_return_cf(u, tau, r, q, params);
        cplx shifted = phi * std::exp(cplx(0.0, 1.0) * u * (x - a));

        double k_pi_over_ba = u;

        double Vk = (2.0 / (b - a)) * strike *
                    (-chi(k_pi_over_ba, a, 0.0, a) + psi(k_idx, k_pi_over_ba, a, 0.0, a));

        double weight = (k_idx == 0) ? 0.5 : 1.0;
        sum += weight * shifted.real() * Vk;
    }

    double price = std::exp(-r * tau) * sum;
    return std::max(price, 0.0);
}

}  // namespace

double bates_price(const EuropeanOption& option, const MarketData& market, const BatesParams& params) {
    double put_price = bates_put_price(market.spot, option.strike, option.time_to_expiry,
                                        market.risk_free_rate, market.dividend_yield, params);

    if (option.type == OptionType::Put) {
        return put_price;
    }

    double disc_spot = market.spot * std::exp(-market.dividend_yield * option.time_to_expiry);
    double disc_strike = option.strike * std::exp(-market.risk_free_rate * option.time_to_expiry);
    return put_price + disc_spot - disc_strike;
}

}  // namespace deriv
