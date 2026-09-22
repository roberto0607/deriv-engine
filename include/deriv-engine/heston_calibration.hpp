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

// Calibrates Heston's 5 parameters (v0, kappa, theta, xi, rho) to a set
// of real market quotes by minimizing pricing error against
// heston_price(). The points can come from a single expiry slice (the
// original, still-common usage -- see tools/calibrate_heston.cpp) or be
// pooled across multiple expiries (see tools/calibrate_heston_multi.cpp):
// nothing in this function's implementation assumes a single expiry --
// each point carries its own time_to_expiry and is priced against it
// individually, so pooling points from several expiries into one call
// simply fits one shared (v0, kappa, theta, xi, rho) against all of them
// jointly. What differs between the two usages is identifiability, not
// code:
//
// Single expiry: kappa (mean-reversion speed) and theta (long-run
// variance) are NOT separately identifiable from one expiry's smile --
// many (kappa, theta) pairs price that one maturity almost identically,
// since what shows up in a single-T price is closer to an integrated
// combination of the two than either one alone. An unregularized fit can
// therefore converge to a mathematically valid but arbitrary-looking
// (kappa, theta) pair while still fitting the smile excellently (see the
// "Heston calibration recovers a synthetic smile..." test, and the
// multi-expiry identifiability test below it, for a direct empirical
// demonstration: single-expiry kappa comes back over 200% off the
// ground truth despite a near-perfect price fit).
//
// Multiple expiries, pooled: because variance mean-reverts at a specific
// (kappa, theta) rate, contracts at different maturities see different
// stages of that reversion -- and fitting all of them at once with one
// shared parameter set pins kappa and theta down uniquely. The
// multi-expiry identifiability test below recovers all 5 parameters
// to machine precision from a deliberately bad initial guess, on
// synthetic (noise-free) data. On the real Deribit chain the picture is
// more honest: constant-parameter Heston cannot fit multiple real
// maturities' smiles as tightly as it fits any one of them alone (see
// tools/calibrate_heston_multi.cpp and docs/numerics.md's "Heston
// calibration: single-expiry vs. multi-expiry" section) -- a genuine,
// well-known limitation of a single-factor stochastic-vol model's term
// structure, not a bug in this calibrator. Resolving that gap properly
// needs a term-structure extension (e.g. piecewise-constant or
// Bergomi-style multi-factor vol), which is out of scope here.
//
// spot and r are passed separately (rather than read off the points)
// because every point in a Deribit snapshot at a given moment shares the
// same underlying_price and risk_free_rate; the caller picks one
// representative point's values (or an average) up front. When pooling
// multiple expiries from the same snapshot, this still holds -- spot and
// r don't vary by expiry within one snapshot.
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
