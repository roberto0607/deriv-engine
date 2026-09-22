#pragma once

#include "heston.hpp"
#include "vol_surface.hpp"
#include <vector>

namespace deriv {

struct HestonCalibrationResult {
    HestonParams params;
    bool converged;
    int iterations;

    // Fit quality, measured two ways:
    //  - rmse_relative_price: what the optimizer actually minimizes (see
    //    calibrate_heston's header comment for why price space, not vol
    //    space).
    //  - rmse_iv_pp: the resulting fit translated into implied-vol terms
    //    (percentage points), computed once after optimization finishes,
    //    for an apples-to-apples comparison with the flat-vol calibration
    //    headline result (see build_vol_surface_from_snapshot).
    double rmse_relative_price;
    double rmse_iv_pp;
    int num_points_used;
};

// Calibrates Heston's 5 parameters (v0, kappa, theta, xi, rho) to a slice
// of real market quotes -- typically all contracts at a single expiry,
// pulled from build_vol_surface_from_snapshot() -- by minimizing pricing
// error against heston_price().
//
// Deliberately scoped to a single expiry slice rather than the whole
// chain at once: Heston's five parameters describe one (v0, kappa, theta,
// xi, rho) instantaneous-variance process, which does not have enough
// free parameters to exactly reproduce every expiry's smile simultaneously
// -- that needs either per-expiry recalibration (what this does, called
// once per slice) or a term-structure extension (e.g. piecewise-constant
// or Bergomi-style multi-factor vol), which is out of scope here. This is
// a deliberate scope decision, not a gap -- see docs/numerics.md.
//
// spot and r are passed separately (rather than read off the points)
// because every point in a single-expiry slice shares the same
// underlying_price and risk_free_rate in a Deribit snapshot; the caller
// picks one representative point's values (or an average) up front.
//
// regularization_weight (default 0.0, i.e. off) adds a small penalty for
// straying from initial_guess in log/atanh-space, alongside the price-fit
// term. Why this exists: kappa and theta are not separately identifiable
// from a single expiry (see the scope note above), so with no
// regularization the optimizer is free to run off along that
// unidentified ridge to whatever (kappa, theta, xi) combination it finds
// first -- fitting the real Deribit chain with weight=0 was observed to
// converge to kappa~37, xi~8, a mathematically valid but implausible-
// looking fit (way outside typical calibrated ranges) purely because
// nothing in the objective discourages it. A small weight (0.001 was
// enough in that case) breaks the tie in favor of solutions that stay
// closer to the initial guess, at a measured cost of a few tenths of a
// percentage point of fit quality -- standard Tikhonov-style
// regularization, not a hack. Left at 0.0 by default so callers that
// want a pure least-squares fit (e.g. the synthetic-recovery test, which
// has no identifiability problem to regularize away) get exactly that.
HestonCalibrationResult calibrate_heston(const std::vector<VolSurfacePoint>& points, double spot, double r,
                                          const HestonParams& initial_guess, double regularization_weight = 0.0);

}  // namespace deriv
