#pragma once

#include "tape2.hpp"
#include "dual.hpp"
#include <cmath>

namespace deriv {

// Second-order sibling of ADouble (adouble.hpp): same operator set,
// same tape-recording behavior, but wraps a Dual instead of a plain
// double so the partials it records onto Tape2 carry forward-mode
// tangents. See tape2.hpp for why this is a separate, isolated type
// rather than a template over ADouble.
class ADouble2 {
public:
    Dual value;
    std::size_t idx;

    // Leaf constructor. `seed_dot` is the forward-mode direction this
    // leaf's own value.dot is seeded with: 1.0 for the ONE input you
    // want the second derivative with respect to (spot, for gamma),
    // 0.0 for every other input (sigma, r, T) -- exactly like seeding a
    // directional derivative in ordinary forward-mode AD, just used
    // here to ride along on top of the reverse-mode pass.
    explicit ADouble2(double v, double seed_dot = 0.0)
        : value(v, seed_dot), idx(get_tape2().push_leaf()) {}

private:
    ADouble2(Dual v, std::size_t existing_idx) : value(v), idx(existing_idx) {}

    friend ADouble2 operator+(const ADouble2& a, const ADouble2& b);
    friend ADouble2 operator-(const ADouble2& a, const ADouble2& b);
    friend ADouble2 operator*(const ADouble2& a, const ADouble2& b);
    friend ADouble2 operator/(const ADouble2& a, const ADouble2& b);
    friend ADouble2 operator+(const ADouble2& a, double k);
    friend ADouble2 operator+(double k, const ADouble2& a);
    friend ADouble2 operator-(const ADouble2& a, double k);
    friend ADouble2 operator-(double k, const ADouble2& a);
    friend ADouble2 operator*(const ADouble2& a, double k);
    friend ADouble2 operator*(double k, const ADouble2& a);
    friend ADouble2 operator/(const ADouble2& a, double k);
    friend ADouble2 operator-(const ADouble2& a);
    friend ADouble2 ad2_exp(const ADouble2& a);
    friend ADouble2 ad2_log(const ADouble2& a);
    friend ADouble2 ad2_sqrt(const ADouble2& a);
    friend ADouble2 ad2_norm_cdf(const ADouble2& a);
};

inline ADouble2 operator+(const ADouble2& a, const ADouble2& b) {
    return ADouble2(a.value + b.value, get_tape2().push_binary(a.idx, Dual(1.0), b.idx, Dual(1.0)));
}
inline ADouble2 operator-(const ADouble2& a, const ADouble2& b) {
    return ADouble2(a.value - b.value, get_tape2().push_binary(a.idx, Dual(1.0), b.idx, Dual(-1.0)));
}
inline ADouble2 operator*(const ADouble2& a, const ADouble2& b) {
    // d(a*b)/da = b, d(a*b)/db = a -- and here b/a are Duals, so these
    // partials carry the forward tangent of b/a along with them; that's
    // the whole mechanism.
    return ADouble2(a.value * b.value, get_tape2().push_binary(a.idx, b.value, b.idx, a.value));
}
inline ADouble2 operator/(const ADouble2& a, const ADouble2& b) {
    Dual result = a.value / b.value;
    Dual d_da = Dual(1.0) / b.value;
    Dual d_db = Dual(-1.0) * a.value / (b.value * b.value);
    return ADouble2(result, get_tape2().push_binary(a.idx, d_da, b.idx, d_db));
}

inline ADouble2 operator+(const ADouble2& a, double k) {
    return ADouble2(a.value + Dual(k), get_tape2().push_unary(a.idx, Dual(1.0)));
}
inline ADouble2 operator+(double k, const ADouble2& a) { return a + k; }

inline ADouble2 operator-(const ADouble2& a, double k) {
    return ADouble2(a.value - Dual(k), get_tape2().push_unary(a.idx, Dual(1.0)));
}
inline ADouble2 operator-(double k, const ADouble2& a) {
    return ADouble2(Dual(k) - a.value, get_tape2().push_unary(a.idx, Dual(-1.0)));
}

inline ADouble2 operator*(const ADouble2& a, double k) {
    return ADouble2(a.value * Dual(k), get_tape2().push_unary(a.idx, Dual(k)));
}
inline ADouble2 operator*(double k, const ADouble2& a) { return a * k; }

inline ADouble2 operator/(const ADouble2& a, double k) {
    return ADouble2(a.value / Dual(k), get_tape2().push_unary(a.idx, Dual(1.0 / k)));
}

inline ADouble2 operator-(const ADouble2& a) {
    return ADouble2(-a.value, get_tape2().push_unary(a.idx, Dual(-1.0)));
}

inline ADouble2 ad2_exp(const ADouble2& a) {
    Dual result = dual_exp(a.value);
    // d(exp(x))/dx = exp(x) -- same as ad_exp, but `result` is now a
    // Dual, so this partial's own .dot is exp(x)*x.dot automatically.
    return ADouble2(result, get_tape2().push_unary(a.idx, result));
}
inline ADouble2 ad2_log(const ADouble2& a) {
    // d(log(x))/dx = 1/x, as a Dual: val = 1/x.val, dot = -x.dot/x.val^2
    // (the forward tangent of that same partial, needed by the second
    // backward pass exactly like every other special function here).
    Dual partial(1.0 / a.value.val, -a.value.dot / (a.value.val * a.value.val));
    Dual result(std::log(a.value.val), a.value.dot / a.value.val);
    return ADouble2(result, get_tape2().push_unary(a.idx, partial));
}

inline ADouble2 ad2_sqrt(const ADouble2& a) {
    Dual result = dual_sqrt(a.value);
    Dual half_over_result = Dual(0.5) / result;
    return ADouble2(result, get_tape2().push_unary(a.idx, half_over_result));
}

// Standard normal CDF, N(x) = 0.5*erfc(-x/sqrt(2)), as a special
// function with a hand-supplied derivative -- the same pattern as
// ad2_exp/ad2_sqrt above (rather than differentiating through erfc's
// own internals, which this codebase doesn't implement from scratch
// anywhere). N'(x) = phi(x), the standard normal PDF, itself
// evaluated on the plain x.value.val the way ad2_exp reuses its own
// result as the reverse-mode partial. This is what lets gamma be
// computed via this second-order tape directly on a smooth closed-
// form pricing formula (Black-Scholes) -- see
// monte_carlo_gamma()/tests/test_main.cpp for why the SAME technique,
// applied instead to a Monte Carlo path's payoff, gives exactly zero
// (a real, documented finding, not a bug -- see docs/phase-log.md
// Phase 11).
inline ADouble2 ad2_norm_cdf(const ADouble2& a) {
    constexpr double inv_sqrt_2pi = 0.3989422804014327;  // 1/sqrt(2*pi)
    double x = a.value.val;
    double cdf_val = 0.5 * std::erfc(-x / std::sqrt(2.0));
    double pdf_val = inv_sqrt_2pi * std::exp(-0.5 * x * x);

    // phi(x) as a Dual: its own value is pdf_val, and its own forward
    // tangent is phi'(x)*x.dot = -x*phi(x)*x.dot (standard normal PDF
    // derivative). This Dual IS the local reverse-mode partial dN/dx.
    Dual phi_as_dual(pdf_val, -x * pdf_val * a.value.dot);
    Dual result(cdf_val, pdf_val * a.value.dot);
    return ADouble2(result, get_tape2().push_unary(a.idx, phi_as_dual));
}

}  // namespace deriv
