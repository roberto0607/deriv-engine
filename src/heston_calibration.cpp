#include "deriv-engine/heston_calibration.hpp"
#include "deriv-engine/implied_vol.hpp"
#include "deriv-engine/market_data.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>

namespace deriv {

namespace {

// Heston's parameters are constrained (v0, kappa, theta, xi > 0;
// -1 < rho < 1), but Nelder-Mead's simplex moves freely through R^5 and
// will happily propose a negative variance or a rho outside [-1, 1] if
// left unconstrained. Rather than clamping (which flattens the objective
// at the boundary and stalls the simplex there), the optimizer works in
// an unconstrained space y and this transform maps y -> a always-valid
// HestonParams: exp() for the four positive parameters, tanh() for rho.
// This is a standard trick for box-constrained Nelder-Mead and keeps
// every point the optimizer ever evaluates physically valid.
HestonParams params_from_unconstrained(const std::array<double, 5>& y) {
    HestonParams p;
    p.v0 = std::exp(y[0]);
    p.kappa = std::exp(y[1]);
    p.theta = std::exp(y[2]);
    p.xi = std::exp(y[3]);
    p.rho = std::tanh(y[4]);
    return p;
}

std::array<double, 5> unconstrained_from_params(const HestonParams& p) {
    return {std::log(p.v0), std::log(p.kappa), std::log(p.theta), std::log(p.xi), std::atanh(std::clamp(p.rho, -0.999, 0.999))};
}

// Generic Nelder-Mead simplex minimizer over R^N. Chosen over a gradient-
// based method (e.g. Levenberg-Marquardt) deliberately: heston_price()
// has no analytic derivative with respect to (v0, kappa, theta, xi, rho)
// implemented anywhere in this codebase (unlike the Greeks, which are
// derivatives with respect to market inputs, not model parameters), and
// numerically differencing a 256-term COS summation 5 times per gradient
// step is both slower and noisier than just not needing a gradient at
// all. Nelder-Mead is the standard choice for exactly this
// (low-dimensional, no gradient, reasonably well-behaved objective)
// situation.
template <std::size_t N>
struct NelderMeadResult {
    std::array<double, N> x;
    double f;
    int iterations;
    bool converged;
};

template <std::size_t N, typename F>
NelderMeadResult<N> nelder_mead(F objective, std::array<double, N> x0, double initial_step = 0.5,
                                 int max_iterations = 2000, double f_tol = 1e-10) {
    // Standard Nelder-Mead coefficients (Nelder & Mead 1965).
    const double alpha = 1.0;  // reflection
    const double gamma = 2.0;  // expansion
    const double rho_c = 0.5;  // contraction
    const double sigma = 0.5;  // shrink

    std::array<std::array<double, N>, N + 1> simplex;
    std::array<double, N + 1> fvals;

    simplex[0] = x0;
    for (std::size_t i = 0; i < N; ++i) {
        simplex[i + 1] = x0;
        simplex[i + 1][i] += initial_step;
    }
    for (std::size_t i = 0; i <= N; ++i) fvals[i] = objective(simplex[i]);

    int iter = 0;
    bool converged = false;

    for (; iter < max_iterations; ++iter) {
        // Sort simplex vertices by objective value, best first.
        std::array<std::size_t, N + 1> order;
        for (std::size_t i = 0; i <= N; ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) { return fvals[a] < fvals[b]; });

        std::array<std::array<double, N>, N + 1> sorted_simplex;
        std::array<double, N + 1> sorted_fvals;
        for (std::size_t i = 0; i <= N; ++i) {
            sorted_simplex[i] = simplex[order[i]];
            sorted_fvals[i] = fvals[order[i]];
        }
        simplex = sorted_simplex;
        fvals = sorted_fvals;

        if (std::abs(fvals[N] - fvals[0]) < f_tol) {
            converged = true;
            break;
        }

        // Centroid of all points except the worst.
        std::array<double, N> centroid{};
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t d = 0; d < N; ++d) centroid[d] += simplex[i][d];
        }
        for (std::size_t d = 0; d < N; ++d) centroid[d] /= static_cast<double>(N);

        auto reflect = [&](double coeff) {
            std::array<double, N> pt;
            for (std::size_t d = 0; d < N; ++d) pt[d] = centroid[d] + coeff * (centroid[d] - simplex[N][d]);
            return pt;
        };

        std::array<double, N> xr = reflect(alpha);
        double fr = objective(xr);

        if (fr < fvals[0]) {
            // Best so far -- try expanding further in the same direction.
            std::array<double, N> xe = reflect(alpha * gamma);
            double fe = objective(xe);
            if (fe < fr) {
                simplex[N] = xe;
                fvals[N] = fe;
            } else {
                simplex[N] = xr;
                fvals[N] = fr;
            }
        } else if (fr < fvals[N - 1]) {
            simplex[N] = xr;
            fvals[N] = fr;
        } else {
            // Reflection didn't help -- contract toward the centroid.
            std::array<double, N> xc;
            for (std::size_t d = 0; d < N; ++d) {
                xc[d] = centroid[d] + rho_c * ((fr < fvals[N] ? xr[d] : simplex[N][d]) - centroid[d]);
            }
            double fc = objective(xc);
            if (fc < std::min(fr, fvals[N])) {
                simplex[N] = xc;
                fvals[N] = fc;
            } else {
                // Contraction also failed -- shrink the whole simplex
                // toward the best point.
                for (std::size_t i = 1; i <= N; ++i) {
                    for (std::size_t d = 0; d < N; ++d) {
                        simplex[i][d] = simplex[0][d] + sigma * (simplex[i][d] - simplex[0][d]);
                    }
                    fvals[i] = objective(simplex[i]);
                }
            }
        }
    }

    // Final sort so x[0]/f[0] is genuinely the best point found.
    std::size_t best = 0;
    for (std::size_t i = 1; i <= N; ++i) {
        if (fvals[i] < fvals[best]) best = i;
    }

    return NelderMeadResult<N>{simplex[best], fvals[best], iter, converged};
}

}  // namespace

HestonCalibrationResult calibrate_heston(const std::vector<VolSurfacePoint>& points, double spot, double r,
                                          const HestonParams& initial_guess, double regularization_weight) {
    std::vector<const VolSurfacePoint*> fit_points;
    for (const auto& p : points) {
        if (p.converged && p.market_price_usd > 0.0) fit_points.push_back(&p);
    }

    HestonCalibrationResult result;
    result.num_points_used = static_cast<int>(fit_points.size());

    if (fit_points.empty()) {
        result.params = initial_guess;
        result.converged = false;
        result.iterations = 0;
        result.rmse_relative_price = 0.0;
        result.rmse_iv_pp = 0.0;
        return result;
    }

    // Objective: root-mean-square *relative* price error, not vol-space
    // error. Relative (not absolute) because Deribit BTC option prices
    // span several orders of magnitude across one expiry's strike range
    // (a deep ITM put and a far OTM call on the same underlying can
    // differ by 1000x in USD price) -- an absolute-error objective would
    // be dominated entirely by the few most expensive (deep ITM)
    // contracts and effectively ignore the OTM wings that carry most of
    // the smile information. Price space (not implied-vol space, which
    // is the more common textbook choice) because it avoids an inverse
    // implied-vol solve inside the objective -- inside a Nelder-Mead loop
    // that already runs a 256-term COS summation per candidate point,
    // that solve would multiply the already-nontrivial calibration cost
    // for comparatively little benefit at this project's scale. The
    // resulting fit is still reported in vol-space terms afterward (see
    // rmse_iv_pp below) for a direct, honest comparison against the
    // flat-vol calibration headline result.
    std::array<double, 5> y0 = unconstrained_from_params(initial_guess);

    auto objective = [&](const std::array<double, 5>& y) -> double {
        HestonParams params = params_from_unconstrained(y);
        MarketData market{spot, r, 0.0, 0.0};

        double sum_sq = 0.0;
        for (const auto* p : fit_points) {
            EuropeanOption option{p->strike, p->time_to_expiry, p->type};
            double model_price = heston_price(option, market, params);
            double rel_err = (model_price - p->market_price_usd) / p->market_price_usd;
            sum_sq += rel_err * rel_err;
        }
        double data_term = sum_sq / static_cast<double>(fit_points.size());

        double reg_term = 0.0;
        if (regularization_weight > 0.0) {
            for (std::size_t i = 0; i < 5; ++i) {
                double diff = y[i] - y0[i];
                reg_term += diff * diff;
            }
            reg_term *= regularization_weight;
        }
        return data_term + reg_term;
    };

    auto opt = nelder_mead<5>(objective, y0, /*initial_step=*/0.5, /*max_iterations=*/3000, /*f_tol=*/1e-12);

    result.params = params_from_unconstrained(opt.x);
    result.converged = opt.converged;
    result.iterations = opt.iterations;
    result.rmse_relative_price = std::sqrt(opt.f);

    // Translate the final fit into implied-vol terms: reprice every
    // fitted contract under the calibrated Heston params, then run it
    // back through this engine's own Newton-Raphson solver (the same one
    // used for the flat-vol headline result) to get a Black-Scholes-
    // equivalent implied vol, and compare that to each contract's
    // independently-solved market implied vol.
    {
        MarketData market{spot, r, 0.0, 0.0};
        double sum_sq_pp = 0.0;
        int n = 0;
        for (const auto* p : fit_points) {
            EuropeanOption option{p->strike, p->time_to_expiry, p->type};
            double model_price = heston_price(option, market, result.params);
            ImpliedVolResult iv = implied_volatility(option, market, model_price);
            if (!iv.converged) continue;
            double diff_pp = (iv.vol - p->solved_iv) * 100.0;
            sum_sq_pp += diff_pp * diff_pp;
            ++n;
        }
        result.rmse_iv_pp = n > 0 ? std::sqrt(sum_sq_pp / n) : 0.0;
    }

    return result;
}

}  // namespace deriv
