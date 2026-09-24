#pragma once

// Shared internals for COS-method (Fang & Oosterlee, 2008) pricers built
// on a characteristic function: currently heston.cpp and bates.cpp. Not
// a public API -- everything here lives in namespace deriv::detail and
// exists purely so bates.cpp can reuse heston.cpp's numerically-hardened
// characteristic function (and both can reuse the same chi/psi quadrature
// coefficients) instead of re-deriving or copy-pasting fragile complex-
// analysis code a second time. See heston.cpp's own comments for the
// full derivation and the specific numerical-stability bugs this
// formulation avoids (catastrophic cancellation in both (A-d) and the
// log((1-g*E)/(1-g)) term).

#include "heston.hpp"
#include <cmath>
#include <complex>

namespace deriv {
namespace detail {

using cplx = std::complex<double>;

// Complex log1p: computes log(1 + z) accurately even when |z| is tiny.
// See heston.cpp's original comment on this function for why the naive
// std::log(1.0 + z) loses precision here.
inline cplx clog1p(cplx z) {
    double zr = z.real();
    double zi = z.imag();
    double re = 0.5 * std::log1p(2.0 * zr + zr * zr + zi * zi);
    double im = std::atan2(zi, 1.0 + zr);
    return cplx(re, im);
}

// Characteristic function of the log-return ln(S_T / S0) under Heston, in
// the "little trap" formulation -- moved here verbatim from heston.cpp so
// bates.cpp can call it with an adjusted drift (r - lambda*k, folding in
// the jump-compensated drift) rather than duplicating this numerically
// careful derivation. See heston.cpp for the full derivation notes.
inline cplx heston_log_return_cf(double u, double tau, double r, double q, const HestonParams& p) {
    const cplx i(0.0, 1.0);
    const cplx iu = i * u;

    cplx rho_xi_iu = p.rho * p.xi * iu;
    cplx A = p.kappa - rho_xi_iu;
    cplx d = std::sqrt(A * A + p.xi * p.xi * (iu + u * u));

    cplx A_minus_d_over_xi_sq = -(iu + u * u) / (A + d);
    cplx A_minus_d = A_minus_d_over_xi_sq * p.xi * p.xi;
    cplx g = A_minus_d / (A + d);

    cplx exp_neg_d_tau = std::exp(-d * tau);

    cplx ratio_minus_1 = g * (1.0 - exp_neg_d_tau) / (1.0 - g);
    cplx log_term = clog1p(ratio_minus_1);

    cplx C = iu * (r - q) * tau + p.kappa * p.theta * tau * A_minus_d_over_xi_sq -
             2.0 * (p.kappa * p.theta / (p.xi * p.xi)) * log_term;

    cplx D = A_minus_d_over_xi_sq * ((1.0 - exp_neg_d_tau) / (1.0 - g * exp_neg_d_tau));

    return std::exp(C + D * p.v0);
}

// chi_k(lower, upper) = integral_{lower}^{upper} e^y cos(k*pi*(y-a)/(b-a)) dy,
// the standard COS-method payoff coefficient (Fang & Oosterlee 2008).
inline double chi(double k_pi_over_ba, double lower, double upper, double a) {
    double denom = 1.0 + k_pi_over_ba * k_pi_over_ba;
    double term1 = std::cos(k_pi_over_ba * (upper - a)) * std::exp(upper) -
                   std::cos(k_pi_over_ba * (lower - a)) * std::exp(lower);
    double term2 = k_pi_over_ba * (std::sin(k_pi_over_ba * (upper - a)) * std::exp(upper) -
                                    std::sin(k_pi_over_ba * (lower - a)) * std::exp(lower));
    return (term1 + term2) / denom;
}

// psi_k(lower, upper) = integral_{lower}^{upper} cos(k*pi*(y-a)/(b-a)) dy.
inline double psi(int k, double k_pi_over_ba, double lower, double upper, double a) {
    if (k == 0) return upper - lower;
    return (std::sin(k_pi_over_ba * (upper - a)) - std::sin(k_pi_over_ba * (lower - a))) / k_pi_over_ba;
}

}  // namespace detail
}  // namespace deriv
