#include "deriv-engine/svi.hpp"
#include <cmath>

namespace deriv {

namespace {

// u = k - m, s = sqrt(u^2 + sigma^2). Exact closed-form first and second
// derivatives of w(k) w.r.t. k, used by both the butterfly condition
// below and (for w' only) nowhere else yet -- kept here rather than in
// the header since callers only ever need w, g, or the two arbitrage
// checks, never the raw derivatives directly.
struct SviDerivatives {
    double w;
    double dw;   // w'(k)
    double d2w;  // w''(k)
};

SviDerivatives svi_derivatives(const SviParams& p, double k) {
    double u = k - p.m;
    double s = std::sqrt(u * u + p.sigma * p.sigma);
    SviDerivatives d;
    d.w = p.a + p.b * (p.rho * u + s);
    d.dw = p.b * (p.rho + u / s);
    d.d2w = p.b * p.sigma * p.sigma / (s * s * s);
    return d;
}

}  // namespace

double svi_total_variance(const SviParams& params, double log_moneyness) {
    double u = log_moneyness - params.m;
    double s = std::sqrt(u * u + params.sigma * params.sigma);
    return params.a + params.b * (params.rho * u + s);
}

double svi_implied_vol(const SviParams& params, double log_moneyness, double time_to_expiry) {
    double w = svi_total_variance(params, log_moneyness);
    if (w <= 0.0 || time_to_expiry <= 0.0) return 0.0;
    return std::sqrt(w / time_to_expiry);
}

double svi_butterfly_condition(const SviParams& params, double log_moneyness) {
    SviDerivatives d = svi_derivatives(params, log_moneyness);
    if (d.w <= 0.0) return -1.0;  // a non-positive total variance is already nonsensical

    double k = log_moneyness;
    double term1 = 1.0 - (k * d.dw) / (2.0 * d.w);
    double term2 = (d.dw / 2.0);
    double g = term1 * term1 - term2 * term2 * (1.0 / d.w + 0.25) + d.d2w / 2.0;
    return g;
}

bool svi_check_butterfly_arbitrage_free(const SviParams& params, double k_min, double k_max, int num_samples,
                                         double tolerance) {
    if (num_samples < 2) return true;
    for (int i = 0; i < num_samples; ++i) {
        double k = k_min + (k_max - k_min) * static_cast<double>(i) / static_cast<double>(num_samples - 1);
        if (svi_butterfly_condition(params, k) < tolerance) return false;
    }
    return true;
}

bool svi_check_calendar_arbitrage_free(const SviParams& params_near, const SviParams& params_far, double k_min,
                                        double k_max, int num_samples, double tolerance) {
    if (num_samples < 2) return true;
    for (int i = 0; i < num_samples; ++i) {
        double k = k_min + (k_max - k_min) * static_cast<double>(i) / static_cast<double>(num_samples - 1);
        double w_near = svi_total_variance(params_near, k);
        double w_far = svi_total_variance(params_far, k);
        if (w_far - w_near < tolerance) return false;
    }
    return true;
}

}  // namespace deriv
