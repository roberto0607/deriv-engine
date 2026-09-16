#pragma once

#include "tape.hpp"
#include <cmath>

namespace deriv {

class ADouble {
public:
    double value;
    std::size_t idx;

    // Constructing from a raw double creates a new LEAF node — use this
    // for genuine inputs (spot, vol, rate) you want derivatives with
    // respect to.
    explicit ADouble(double v) : value(v), idx(get_tape().push_leaf()) {}

private:
    // Internal constructor used by operators to wrap an already-computed
    // value with an already-pushed tape index, without creating a
    // redundant new leaf.
    ADouble(double v, std::size_t existing_idx) : value(v), idx(existing_idx) {}

    friend ADouble operator+(const ADouble& a, const ADouble& b);
    friend ADouble operator-(const ADouble& a, const ADouble& b);
    friend ADouble operator*(const ADouble& a, const ADouble& b);
    friend ADouble operator/(const ADouble& a, const ADouble& b);
    friend ADouble operator+(const ADouble& a, double k);
    friend ADouble operator+(double k, const ADouble& a);
    friend ADouble operator-(const ADouble& a, double k);
    friend ADouble operator-(double k, const ADouble& a);
    friend ADouble operator*(const ADouble& a, double k);
    friend ADouble operator*(double k, const ADouble& a);
    friend ADouble operator/(const ADouble& a, double k);
    friend ADouble operator/(double k, const ADouble& a);
    friend ADouble operator-(const ADouble& a);
    friend ADouble ad_exp(const ADouble& a);
    friend ADouble ad_log(const ADouble& a);
    friend ADouble ad_sqrt(const ADouble& a);
};

// ADouble op ADouble — both operands are recorded nodes on the tape
inline ADouble operator+(const ADouble& a, const ADouble& b) {
    return ADouble(a.value + b.value, get_tape().push_binary(a.idx, 1.0, b.idx, 1.0));
}
inline ADouble operator-(const ADouble& a, const ADouble& b) {
    return ADouble(a.value - b.value, get_tape().push_binary(a.idx, 1.0, b.idx, -1.0));
}
inline ADouble operator*(const ADouble& a, const ADouble& b) {
    return ADouble(a.value * b.value, get_tape().push_binary(a.idx, b.value, b.idx, a.value));
}
inline ADouble operator/(const ADouble& a, const ADouble& b) {
    double result = a.value / b.value;
    // d(a/b)/da = 1/b,  d(a/b)/db = -a/b^2
    return ADouble(result, get_tape().push_binary(a.idx, 1.0 / b.value, b.idx, -a.value / (b.value * b.value)));
}

// ADouble op double — the double is a plain CONSTANT, not a new leaf.
// We fold it directly into the local derivative instead of taping it,
// since we never want a derivative with respect to a literal constant.
inline ADouble operator+(const ADouble& a, double k) {
    return ADouble(a.value + k, get_tape().push_unary(a.idx, 1.0));
}
inline ADouble operator+(double k, const ADouble& a) { return a + k; }

inline ADouble operator-(const ADouble& a, double k) {
    return ADouble(a.value - k, get_tape().push_unary(a.idx, 1.0));
}
inline ADouble operator-(double k, const ADouble& a) {
    return ADouble(k - a.value, get_tape().push_unary(a.idx, -1.0));
}

inline ADouble operator*(const ADouble& a, double k) {
    return ADouble(a.value * k, get_tape().push_unary(a.idx, k));
}
inline ADouble operator*(double k, const ADouble& a) { return a * k; }

inline ADouble operator/(const ADouble& a, double k) {
    return ADouble(a.value / k, get_tape().push_unary(a.idx, 1.0 / k));
}
inline ADouble operator/(double k, const ADouble& a) {
    double result = k / a.value;
    return ADouble(result, get_tape().push_unary(a.idx, -k / (a.value * a.value)));
}

inline ADouble operator-(const ADouble& a) {
    return ADouble(-a.value, get_tape().push_unary(a.idx, -1.0));
}

// Named ad_exp/ad_log/ad_sqrt (not exp/log/sqrt) to avoid silently
// shadowing std::exp/std::log/std::sqrt for plain doubles elsewhere
// in the codebase.
inline ADouble ad_exp(const ADouble& a) {
    double result = std::exp(a.value);
    return ADouble(result, get_tape().push_unary(a.idx, result));  // d(exp(x))/dx = exp(x)
}
inline ADouble ad_log(const ADouble& a) {
    return ADouble(std::log(a.value), get_tape().push_unary(a.idx, 1.0 / a.value));
}
inline ADouble ad_sqrt(const ADouble& a) {
    double result = std::sqrt(a.value);
    return ADouble(result, get_tape().push_unary(a.idx, 0.5 / result));  // d(sqrt(x))/dx = 1/(2*sqrt(x))
}

}  // namespace deriv