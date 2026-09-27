#include "deriv-engine/svi_calibration.hpp"
#include "deriv-engine/optimizers.hpp"
#include <algorithm>
#include <array>
#include <cmath>

namespace deriv {

namespace {

// SVI's parameters are constrained (b >= 0, sigma > 0, -1 < rho < 1; a is
// technically free but see the header note on why this doesn't need its
// own bound). Same box-constrained Nelder-Mead trick heston_calibration.cpp
// uses: optimize over an unconstrained y and map it through exp()/tanh()
// so every candidate the simplex ever proposes is automatically valid.
SviParams params_from_unconstrained(const std::array<double, 5>& y) {
    SviParams p;
    p.a = y[0];
    p.b = std::exp(y[1]);
    p.rho = std::tanh(y[2]);
    p.m = y[3];
    p.sigma = std::exp(y[4]);
    return p;
}

std::array<double, 5> unconstrained_from_params(const SviParams& p) {
    return {p.a, std::log(std::max(p.b, 1e-6)), std::atanh(std::clamp(p.rho, -0.999, 0.999)), p.m,
            std::log(std::max(p.sigma, 1e-6))};
}

}  // namespace

SviCalibrationResult calibrate_svi(const std::vector<VolSurfacePoint>& points_for_one_expiry, double spot, double r,
                                    double time_to_expiry) {
    double forward = spot * std::exp(r * time_to_expiry);  // q = 0 -- see instrument.hpp's BTC-scope note

    struct FitPoint {
        double log_moneyness;
        double w_market;  // (solved_iv)^2 * T -- the quantity SVI actually fits, see header comment
    };
    std::vector<FitPoint> fit_points;
    for (const auto& p : points_for_one_expiry) {
        if (!p.converged || p.solved_iv <= 0.0) continue;
        double k = std::log(p.strike / forward);
        double w = p.solved_iv * p.solved_iv * time_to_expiry;
        fit_points.push_back({k, w});
    }

    SviCalibrationResult result;
    result.forward = forward;
    result.time_to_expiry = time_to_expiry;
    result.num_points_used = static_cast<int>(fit_points.size());

    if (fit_points.size() < 5) {
        // Fewer points than free parameters -- refuse to fit rather than
        // return an overfit, meaningless curve. 5 real strikes at one
        // expiry is a low bar; a genuinely thin expiry on the live chain
        // should report "not enough data," not a curve nobody should
        // trust.
        result.params = SviParams{0.0, 0.0, 0.0, 0.0, 0.1};
        result.converged = false;
        result.iterations = 0;
        result.rmse_iv_pp = 0.0;
        result.butterfly_arbitrage_free = false;
        return result;
    }

    // Initial guess: read a0 (the level) straight off the point closest
    // to at-the-money, and start with a mild negative skew (rho0 = -0.3,
    // the common starting point in the SVI literature for equity-like
    // smiles) and a modest curvature -- Nelder-Mead only needs a
    // reasonable starting region, not a good guess, since it's a
    // low-dimensional, well-behaved-in-practice objective.
    auto atm_it = std::min_element(fit_points.begin(), fit_points.end(), [](const FitPoint& x, const FitPoint& y) {
        return std::abs(x.log_moneyness) < std::abs(y.log_moneyness);
    });
    double a0 = atm_it->w_market;
    SviParams initial_guess{a0, 0.1, -0.3, 0.0, 0.1};
    std::array<double, 5> y0 = unconstrained_from_params(initial_guess);

    auto objective = [&](const std::array<double, 5>& y) -> double {
        SviParams params = params_from_unconstrained(y);
        double sum_sq = 0.0;
        for (const auto& fp : fit_points) {
            double w_model = svi_total_variance(params, fp.log_moneyness);
            double diff = w_model - fp.w_market;
            sum_sq += diff * diff;
        }
        return sum_sq / static_cast<double>(fit_points.size());
    };

    auto opt = nelder_mead<5>(objective, y0, /*initial_step=*/0.3, /*max_iterations=*/5000, /*f_tol=*/1e-12);

    result.params = params_from_unconstrained(opt.x);
    result.converged = opt.converged;
    result.iterations = opt.iterations;

    double sum_sq_pp = 0.0;
    int n = 0;
    for (const auto& fp : fit_points) {
        double vol_model = svi_implied_vol(result.params, fp.log_moneyness, time_to_expiry);
        double vol_market = std::sqrt(fp.w_market / time_to_expiry);
        if (vol_model <= 0.0) continue;
        double diff_pp = (vol_model - vol_market) * 100.0;
        sum_sq_pp += diff_pp * diff_pp;
        ++n;
    }
    result.rmse_iv_pp = n > 0 ? std::sqrt(sum_sq_pp / n) : 0.0;

    // Sample the fitted curve well beyond the real strikes it was fit to
    // ([-2, 2] in log-moneyness is roughly a 7x-move range either way --
    // real BTC strikes rarely span that far) so the check catches
    // wing behavior the fit itself never saw, not just interpolation
    // between real points.
    result.butterfly_arbitrage_free = svi_check_butterfly_arbitrage_free(result.params);

    return result;
}

}  // namespace deriv
