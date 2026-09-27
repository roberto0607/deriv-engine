#pragma once

#include <vector>

namespace deriv {

// Raw SVI ("Stochastic Volatility Inspired", Gatheral 2004) parameters for
// one expiry's implied-total-variance smile, as a function of
// log-moneyness k = log(K / F) (K = strike, F = forward price).
//
// This is a curve-fitting parameterization, not a stochastic-process
// model like Heston/Bates: it doesn't describe how the underlying moves,
// it directly describes the shape of the smile itself, fit independently
// at each expiry to whatever the real market is quoting right now. This
// is deliberately a different philosophy from Heston (Phase 9) and Bates
// (Phase 14) already in this engine -- see svi_calibration.hpp for why
// both are worth having rather than one replacing the other.
struct SviParams {
    double a;      // overall variance level
    double b;      // angle between the left/right wings (b >= 0)
    double rho;    // wing asymmetry / skew, -1 < rho < 1
    double m;      // horizontal shift of the smile's minimum
    double sigma;  // curvature at the minimum (sigma > 0)
};

// Gatheral's raw SVI formula for total implied variance (NOT vol; total
// variance = vol^2 * T) at log-moneyness k:
//
//   w(k) = a + b * ( rho * (k - m) + sqrt((k - m)^2 + sigma^2) )
//
// Total variance (rather than vol itself) is SVI's natural unit -- it's
// what makes the calendar no-arbitrage check below a simple inequality
// instead of a T-dependent rescaling, and it's what the calibration
// objective in svi_calibration.cpp fits directly.
double svi_total_variance(const SviParams& params, double log_moneyness);

// Converts svi_total_variance's output at a given expiry T (in years)
// back into an ordinary annualized volatility: vol = sqrt(w(k) / T).
// This is the number every other card on this page already understands
// (Black-Scholes, Heston, Bates all take a vol, not a total variance).
double svi_implied_vol(const SviParams& params, double log_moneyness, double time_to_expiry);

// Gatheral's butterfly (static) no-arbitrage condition: a fitted SVI
// curve is butterfly-arbitrage-free at k iff
//
//   g(k) = (1 - k*w'(k)/(2*w(k)))^2 - (w'(k)/2)^2 * (1/w(k) + 1/4) + w''(k)/2 >= 0
//
// A negative g(k) at some strike means the fitted curve implies a
// negative risk-neutral density there -- i.e. a mathematically
// nonsensical smile, not just an inaccurate one (a real butterfly spread
// built from three of your own model's option prices would show a
// negative cost, which is a genuine arbitrage against your own model).
// w'/w'' are computed by exact closed-form differentiation of
// svi_total_variance (not finite differences), since SVI's derivatives
// have a simple closed form.
double svi_butterfly_condition(const SviParams& params, double log_moneyness);

// Checks butterfly_condition >= 0 (with a small numerical tolerance) at
// `num_samples` points spanning [k_min, k_max], and returns whether every
// sampled point passed. This is a sampled check, not a proof for every
// real k in the interval -- see svi_calibration.hpp's header comment for
// why that's the honest, standard way this check is done in practice.
bool svi_check_butterfly_arbitrage_free(const SviParams& params, double k_min = -2.0, double k_max = 2.0,
                                         int num_samples = 200, double tolerance = -1e-9);

// Gatheral & Jacquier's sufficient calendar-spread no-arbitrage condition
// between two expiries T1 < T2 sharing the same forward: total variance
// must not decrease with time at any shared log-moneyness, i.e.
// w(k; params_far) >= w(k; params_near) for every k. If it does decrease
// somewhere, you could in principle sell the far-expiry option and buy
// the near-expiry option at that strike and lock in a model-implied
// profit regardless of what BTC does -- a real arbitrage against the
// fitted surface, not just two independently-plausible-looking curves
// that happen to disagree.
bool svi_check_calendar_arbitrage_free(const SviParams& params_near, const SviParams& params_far,
                                        double k_min = -2.0, double k_max = 2.0, int num_samples = 200,
                                        double tolerance = -1e-9);

}  // namespace deriv
