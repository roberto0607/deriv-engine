#pragma once

#include <algorithm>
#include <array>
#include <cstddef>

namespace deriv {

// Generic Nelder-Mead simplex minimizer over R^N. Originally written
// inside heston_calibration.cpp for Heston's 5-parameter fit; pulled out
// here (Phase 21) so SVI's own 5-parameter fit (svi_calibration.cpp)
// reuses the exact same optimizer instead of a second, subtly different
// copy. Nothing about this implementation is specific to Heston or SVI --
// it minimizes any objective function over R^N with no derivative
// required, which is why it works for both: neither heston_price() (a
// 256-term COS summation) nor svi_total_variance() (closed-form, but the
// calibration objective still isn't smooth enough to trust a gradient
// through a plain finite-difference approximation) has an analytic
// derivative with respect to its own model parameters implemented
// anywhere in this codebase.
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

}  // namespace deriv
