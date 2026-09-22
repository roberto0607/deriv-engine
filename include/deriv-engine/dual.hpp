#pragma once

#include <cmath>

namespace deriv {

// A forward-mode dual number: val + dot*eps, with eps^2 = 0. Standard
// forward-mode AD -- propagating (val, dot) through arithmetic via the
// product/quotient rule gives dot = d(val)/d(seed direction) exactly,
// for whatever direction the seed was set to (see push_leaf below).
//
// This exists so Tape2/ADouble2 (adouble2.hpp) can do "forward-over-
// reverse" AD: running the EXISTING reverse-mode tape machinery
// (tape.hpp/adouble.hpp), but with every value and every local partial
// derivative replaced by a Dual instead of a plain double. Reverse mode
// alone gives first derivatives (the Greeks already computed by
// monte_carlo_greeks); layering forward-mode underneath it via Dual
// gives the derivative *of those first derivatives*, i.e. gamma =
// d(delta)/d(spot), without literally differentiating the backward
// pass by hand (a much larger undertaking -- see the note in
// tape2.hpp on why this is the standard technique for exactly this
// situation).
struct Dual {
    double val = 0.0;
    double dot = 0.0;

    Dual() = default;
    // Implicit from a plain double: a genuine CONSTANT, dot = 0 (its
    // value doesn't change as we move in the seeded direction).
    Dual(double v) : val(v), dot(0.0) {}
    Dual(double v, double d) : val(v), dot(d) {}
};

inline Dual operator+(Dual a, Dual b) { return Dual(a.val + b.val, a.dot + b.dot); }
inline Dual operator-(Dual a, Dual b) { return Dual(a.val - b.val, a.dot - b.dot); }
inline Dual operator-(Dual a) { return Dual(-a.val, -a.dot); }

// Product rule: d(a*b) = a*db + b*da.
inline Dual operator*(Dual a, Dual b) { return Dual(a.val * b.val, a.dot * b.val + a.val * b.dot); }

// Quotient rule: d(a/b) = (da*b - a*db) / b^2.
inline Dual operator/(Dual a, Dual b) {
    double val = a.val / b.val;
    double dot = (a.dot * b.val - a.val * b.dot) / (b.val * b.val);
    return Dual(val, dot);
}

// d(exp(u))/deps = exp(u) * du/deps -- chain rule, mirrors ad_exp in adouble.hpp.
inline Dual dual_exp(Dual a) {
    double e = std::exp(a.val);
    return Dual(e, e * a.dot);
}

// d(sqrt(u))/deps = du/deps / (2*sqrt(u)).
inline Dual dual_sqrt(Dual a) {
    double s = std::sqrt(a.val);
    return Dual(s, a.dot / (2.0 * s));
}

}  // namespace deriv
