#pragma once

#include "svi.hpp"
#include "vol_surface.hpp"
#include <vector>

namespace deriv {

struct SviCalibrationResult {
    SviParams params;
    bool converged;
    int iterations;

    double forward;               // F = spot * exp(r * T), q = 0 (see instrument.hpp's BTC-scope note)
    double time_to_expiry;        // the single expiry this fit covers
    double rmse_iv_pp;            // fit quality in implied-vol percentage points
    int num_points_used;
    bool butterfly_arbitrage_free;  // sampled check, see svi.hpp
};

// Fits one expiry's worth of real market points to a single SVI curve
// by minimizing squared error in TOTAL VARIANCE space (w = vol^2 * T),
// not vol space and not price space. Three reasons total-variance space
// specifically, not either of this engine's other two calibration
// choices:
//   - Not price space (heston_calibration.cpp's choice): SVI has no
//     pricing formula of its own -- it parameterizes a smile, not a
//     process, so there's no model price to compare against a market
//     price in the first place. The only quantity SVI's formula produces
//     is total variance.
//   - Not vol space directly: vol = sqrt(w/T), so minimizing squared vol
//     error and minimizing squared total-variance error are NOT the same
//     objective (the sqrt is nonlinear) -- fitting in w-space is what
//     makes the butterfly/calendar arbitrage formulas in svi.hpp exact
//     and is the standard convention in the SVI literature (Gatheral &
//     Jacquier), so this calibrator matches that convention rather than
//     inventing its own.
//
// points must all share the same expiry (same time_to_expiry, up to
// floating-point noise from real timestamp parsing) -- this fits ONE
// SVI curve per call, by design: real Deribit expiries genuinely have
// different-looking smiles, and forcing one shared (a,b,rho,m,sigma)
// across all of them the way Heston's multi-expiry calibration
// deliberately does (see heston_calibration.hpp) would defeat SVI's
// entire purpose, which is to match each expiry's real smile shape
// exactly rather than explain multiple expiries from one
// stochastic-process story. Call this once per expiry present in a live
// chain to build a full surface (see docs/index.html's SVI section for
// how the live demo does exactly that).
//
// spot and r are passed the same way heston_calibration.hpp's
// calibrate_heston takes them: one representative value per snapshot,
// since every point at a given moment shares them.
SviCalibrationResult calibrate_svi(const std::vector<VolSurfacePoint>& points_for_one_expiry, double spot, double r,
                                    double time_to_expiry);

}  // namespace deriv
