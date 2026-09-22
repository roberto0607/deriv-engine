#include "deriv-engine/heston.hpp"
#include <algorithm>
#include <cmath>
#include <complex>

namespace deriv {

namespace {

using cplx = std::complex<double>;

// Complex log1p: computes log(1 + z) accurately even when |z| is tiny,
// the way std::log1p does for reals. std::log(1.0 + z) loses precision
// here because "1.0 + z" rounds to the nearest double *before* the log
// is taken -- if z is, say, O(1e-12), most or all of its digits are
// gone before log() ever sees them. This splits log(1+z) into a real
// part computed via the real std::log1p on |1+z|^2-1 (itself expanded
// so it doesn't re-introduce the same cancellation: |1+z|^2 - 1 =
// 2*Re(z) + |z|^2, which needs no subtraction of near-equal numbers)
// and an imaginary part via atan2, which is well-conditioned on its
// own. See the note on heston_log_return_cf() below for why this
// specific cancellation matters here.
cplx clog1p(cplx z) {
    double zr = z.real();
    double zi = z.imag();
    double re = 0.5 * std::log1p(2.0 * zr + zr * zr + zi * zi);
    double im = std::atan2(zi, 1.0 + zr);
    return cplx(re, im);
}

// Characteristic function of the log-return ln(S_T / S0) under Heston,
// in the "little trap" formulation. This is the piece that's genuinely
// easy to get subtly wrong (a sign flip in g(u), or picking the wrong
// branch of d(u)), so its correctness is verified two independent
// ways elsewhere in this codebase: (1) the deterministic-variance
// limit must collapse exactly to Black-Scholes (see
// tests/test_main.cpp), and (2) an independent Monte Carlo simulation
// of the actual Heston SDEs must agree with this closed-form price.
cplx heston_log_return_cf(double u, double tau, double r, double q, const HestonParams& p) {
    const cplx i(0.0, 1.0);
    const cplx iu = i * u;

    cplx rho_xi_iu = p.rho * p.xi * iu;
    cplx A = p.kappa - rho_xi_iu;
    cplx d = std::sqrt(A * A + p.xi * p.xi * (iu + u * u));

    // (A - d) is computed here as -(iu + u^2) / (A + d), NOT as a
    // direct subtraction of A and d. This matters: A and d are both
    // O(kappa) and nearly equal whenever xi is small-to-moderate (not
    // just in the xi->0 limit), so computing (A - d) as a literal
    // difference loses precision to catastrophic cancellation --
    // confirmed empirically (this cost ~0.5-6% pricing error at
    // xi=0.01-0.5 before this fix, verified against the exact
    // deterministic-variance Black-Scholes limit). The identity
    // A^2 - d^2 = -xi^2 (iu + u^2) lets (A-d) be computed as
    // (A^2-d^2)/(A+d) instead -- a sum in the denominator, not a
    // difference -- which stays numerically stable across the whole
    // parameter range. It also cleanly cancels the xi^2 that
    // multiplies (A-d) everywhere it's used below.
    cplx A_minus_d_over_xi_sq = -(iu + u * u) / (A + d);
    cplx A_minus_d = A_minus_d_over_xi_sq * p.xi * p.xi;
    cplx g = A_minus_d / (A + d);

    cplx exp_neg_d_tau = std::exp(-d * tau);

    // The log((1 - g*E)/(1 - g)) term below (E = exp(-d*tau)) is a
    // second, independent cancellation hazard from the (A-d) one
    // above, and it survives even after that fix: g = O(xi^2) is tiny
    // whenever xi is small-to-moderate (not just at xi->0), which
    // makes the ratio (1-g*E)/(1-g) extremely close to 1. Evaluating
    // std::log() of a value that's within ~1e-10 to 1e-15 of 1.0
    // rounds away most or all of the meaningful digits *before* log()
    // runs, since "1.0 - g*E" and "1.0 - g" both round to (or near) 1.0
    // in double precision.
    //
    // This CF error is normally invisible -- it's a tiny relative error
    // in an intermediate quantity that gets lost in the noise for
    // ordinary (near-the-money) options. It stopped being invisible
    // for deep out-of-the-money/in-the-money options priced through
    // put-call parity (see heston_put_price() below): parity derives a
    // *small* target price as the difference of two *large* numbers
    // (the put and the discounted forward), so the target's accuracy
    // is bounded by the *absolute* accuracy of the put, not its
    // relative accuracy. Confirmed empirically: for
    // spot=100000/strike=300000/90d/60%vol (strike = 3x spot, realistic
    // for a crypto options chain), this cancellation alone accounted
    // for an 82% pricing error in the deterministic-variance limit
    // (which must equal Black-Scholes exactly) -- fixed to ~0.001% by
    // rewriting the ratio as 1 + g*(1-E)/(1-g) and evaluating its log
    // via clog1p() above, which needs no subtraction of near-equal
    // values at all.
    cplx ratio_minus_1 = g * (1.0 - exp_neg_d_tau) / (1.0 - g);
    cplx log_term = clog1p(ratio_minus_1);

    cplx C = iu * (r - q) * tau + p.kappa * p.theta * tau * A_minus_d_over_xi_sq -
             2.0 * (p.kappa * p.theta / (p.xi * p.xi)) * log_term;

    cplx D = A_minus_d_over_xi_sq * ((1.0 - exp_neg_d_tau) / (1.0 - g * exp_neg_d_tau));

    return std::exp(C + D * p.v0);
}

// chi_k(lower, upper) = integral_{lower}^{upper} e^y cos(k*pi*(y-a)/(b-a)) dy,
// the standard COS-method payoff coefficient (Fang & Oosterlee 2008).
double chi(double k_pi_over_ba, double lower, double upper, double a) {
    double denom = 1.0 + k_pi_over_ba * k_pi_over_ba;
    double term1 = std::cos(k_pi_over_ba * (upper - a)) * std::exp(upper) -
                   std::cos(k_pi_over_ba * (lower - a)) * std::exp(lower);
    double term2 = k_pi_over_ba * (std::sin(k_pi_over_ba * (upper - a)) * std::exp(upper) -
                                    std::sin(k_pi_over_ba * (lower - a)) * std::exp(lower));
    return (term1 + term2) / denom;
}

// psi_k(lower, upper) = integral_{lower}^{upper} cos(k*pi*(y-a)/(b-a)) dy.
double psi(int k, double k_pi_over_ba, double lower, double upper, double a) {
    if (k == 0) return upper - lower;
    return (std::sin(k_pi_over_ba * (upper - a)) - std::sin(k_pi_over_ba * (lower - a))) / k_pi_over_ba;
}

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
