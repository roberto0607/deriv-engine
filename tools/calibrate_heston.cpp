// Calibrates Heston to a single expiry slice of the committed real Deribit
// snapshot and emits the result as JSON on stdout -- the same pattern as
// gen_stats.cpp, and for the same reason: never hand-type these numbers
// into docs/index.html. This program (run by CI, see pages.yml) is the
// only source of truth for the calibrated (kappa, xi, rho) the live demo's
// Heston card uses in place of illustrative placeholder values.
//
// Only kappa/xi/rho are consumed by the demo -- v0/theta stay derived from
// the calculator's own flat-vol input (vol^2), same as before calibration,
// since those track "today's ATM vol" rather than a structural model
// parameter. See heston_calibration.hpp for why a single expiry slice
// (not the whole chain) is what gets calibrated, and why moneyness
// filtering is applied before fitting.

#include "deriv-engine/heston_calibration.hpp"
#include "deriv-engine/vol_surface.hpp"
#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>

int main() {
    auto all_points = deriv::build_vol_surface_from_snapshot("data/btc_chain_snapshot.json", 2026, 9, 17);

    const double target_T = 90.0 / 365.0;
    double chosen_key = -1.0;
    {
        std::vector<std::pair<double, int>> expiry_counts;
        for (const auto& p : all_points) {
            double key = std::round(p.time_to_expiry * 3650.0) / 3650.0;
            bool found = false;
            for (auto& kv : expiry_counts) {
                if (std::abs(kv.first - key) < 1e-9) {
                    kv.second++;
                    found = true;
                    break;
                }
            }
            if (!found) expiry_counts.push_back({key, 1});
        }
        double best_dist = 1e18;
        for (const auto& kv : expiry_counts) {
            if (kv.second < 10) continue;
            double dist = std::abs(kv.first - target_T);
            if (dist < best_dist) {
                best_dist = dist;
                chosen_key = kv.first;
            }
        }
    }

    if (chosen_key <= 0.0) {
        std::fprintf(stderr, "calibrate_heston: no expiry with enough contracts found\n");
        return 1;
    }

    std::vector<deriv::VolSurfacePoint> slice;
    double spot = 0.0, r = 0.0;
    for (const auto& p : all_points) {
        double key = std::round(p.time_to_expiry * 3650.0) / 3650.0;
        if (std::abs(key - chosen_key) < 1e-9 && p.converged && p.moneyness > 0.7 && p.moneyness < 1.4) {
            slice.push_back(p);
            spot = p.underlying_price;
            r = p.risk_free_rate;
        }
    }

    if (slice.size() < 10) {
        std::fprintf(stderr, "calibrate_heston: too few contracts in the filtered slice (%zu)\n", slice.size());
        return 1;
    }

    double atm_iv = 0.6;
    double best_moneyness_dist = 1e18;
    for (const auto& p : slice) {
        double dist = std::abs(p.moneyness - 1.0);
        if (dist < best_moneyness_dist) {
            best_moneyness_dist = dist;
            atm_iv = p.solved_iv;
        }
    }

    deriv::HestonParams guess{atm_iv * atm_iv, 2.0, atm_iv * atm_iv, 0.7, -0.5};

    // Regularized (weight=0.001): kappa and theta aren't separately
    // identifiable from one expiry (see heston_calibration.hpp), and an
    // unregularized fit to this real chain was observed to converge to
    // kappa~37/xi~8 -- a valid fit, but nowhere near a typical calibrated
    // Heston parameter range, and not something worth putting on a public
    // demo unqualified. Regularizing toward the (reasonable, ATM-vol-
    // based) initial guess costs a few tenths of a percentage point of
    // fit quality in exchange for kappa/xi landing in a plausible range.
    deriv::HestonCalibrationResult result = deriv::calibrate_heston(slice, spot, r, guess, /*regularization_weight=*/0.001);

    std::printf("{\n");
    std::printf("  \"expiry_days\": %.1f,\n", chosen_key * 365.0);
    std::printf("  \"spot\": %.2f,\n", spot);
    std::printf("  \"contracts_used\": %d,\n", result.num_points_used);
    std::printf("  \"converged\": %s,\n", result.converged ? "true" : "false");
    std::printf("  \"v0\": %.6f,\n", result.params.v0);
    std::printf("  \"kappa\": %.6f,\n", result.params.kappa);
    std::printf("  \"theta\": %.6f,\n", result.params.theta);
    std::printf("  \"xi\": %.6f,\n", result.params.xi);
    std::printf("  \"rho\": %.6f,\n", result.params.rho);
    std::printf("  \"rmse_iv_pp\": %.4f,\n", result.rmse_iv_pp);
    std::printf("  \"rmse_relative_price\": %.6f\n", result.rmse_relative_price);
    std::printf("}\n");

    return 0;
}
