// Calibrates Heston jointly across multiple real expiries from the
// committed Deribit snapshot, instead of the single-expiry slice
// tools/calibrate_heston.cpp uses -- see heston_calibration.hpp's header
// comment for why pooling expiries resolves the kappa/theta
// non-identifiability that a single expiry's smile suffers from.
//
// This program is a diagnostic / research tool, not (yet) wired into the
// live demo's Heston card the way calibrate_heston.cpp is: the finding it
// produces is a real one, but not an unambiguous improvement to headline
// on a public demo -- see the comment above the fit-quality output below.
//
// Run from the repo root: build/calibrate_heston_multi

#include "deriv-engine/heston_calibration.hpp"
#include "deriv-engine/vol_surface.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <utility>
#include <vector>

namespace {

// Picks the expiry (rounded to ~0.1 day) nearest each target day count
// that has at least min_contracts contracts in [0.7, 1.4] moneyness --
// same liquidity/moneyness filter tools/calibrate_heston.cpp uses, for
// the same reason (thin far-wing quotes carry more noise than signal
// under a relative-price objective).
std::vector<double> choose_expiry_keys(const std::vector<deriv::VolSurfacePoint>& all_points,
                                        const std::vector<double>& target_days, int min_contracts) {
    std::vector<std::pair<double, int>> expiry_counts;
    for (const auto& p : all_points) {
        if (!p.converged || p.market_price_usd <= 0.0 || p.moneyness < 0.7 || p.moneyness > 1.4) continue;
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

    std::vector<double> chosen;
    for (double target_day : target_days) {
        double target_T = target_day / 365.0;
        double best_dist = 1e18, best_key = -1.0;
        for (const auto& kv : expiry_counts) {
            if (kv.second < min_contracts) continue;
            double dist = std::abs(kv.first - target_T);
            if (dist < best_dist) {
                best_dist = dist;
                best_key = kv.first;
            }
        }
        if (best_key > 0.0) chosen.push_back(best_key);
    }
    return chosen;
}

}  // namespace

int main() {
    auto all_points = deriv::build_vol_surface_from_snapshot("data/btc_chain_snapshot.json", 2026, 9, 17);

    // Two nearby expiries (~43d, ~99d): close enough in maturity that a
    // constant-parameter Heston fit still has a fighting chance, while
    // still far enough apart to carry real term-structure information.
    // (Adding farther-dated expiries -- 190d, 281d -- was tried during
    // development and pools more identifiability signal, but the fit
    // quality degrades further as a result; see docs/numerics.md.)
    auto keys = choose_expiry_keys(all_points, {43.0, 99.0}, /*min_contracts=*/10);
    if (keys.size() < 2) {
        std::fprintf(stderr, "calibrate_heston_multi: could not find 2 expiries with enough contracts\n");
        return 1;
    }

    std::vector<deriv::VolSurfacePoint> pooled;
    double spot = 0.0, r = 0.0;
    std::vector<std::pair<double, int>> expiry_report;
    for (double key : keys) {
        int count = 0;
        for (const auto& p : all_points) {
            double k = std::round(p.time_to_expiry * 3650.0) / 3650.0;
            if (std::abs(k - key) < 1e-9 && p.converged && p.moneyness > 0.7 && p.moneyness < 1.4) {
                pooled.push_back(p);
                spot = p.underlying_price;
                r = p.risk_free_rate;
                ++count;
            }
        }
        expiry_report.push_back({key * 365.0, count});
    }

    if (pooled.size() < 20) {
        std::fprintf(stderr, "calibrate_heston_multi: too few pooled contracts (%zu)\n", pooled.size());
        return 1;
    }

    // ATM iv from the shortest-dated slice as the vol-level guess (same
    // pattern as tools/calibrate_heston.cpp).
    double atm_iv = 0.6, best_dist = 1e18;
    for (const auto& p : pooled) {
        if (std::abs(p.time_to_expiry - keys[0]) > 1e-9) continue;
        double dist = std::abs(p.moneyness - 1.0);
        if (dist < best_dist) {
            best_dist = dist;
            atm_iv = p.solved_iv;
        }
    }
    deriv::HestonParams guess{atm_iv * atm_iv, 2.0, atm_iv * atm_iv, 0.7, -0.5};

    // Regularized (weight=0.001), matching tools/calibrate_heston.cpp:
    // even with multiple expiries pooled, an unregularized fit on this
    // real (noisy, imperfectly-Heston) chain was observed during
    // development to run off toward v0 -> 0 with an implausible kappa/xi
    // -- the pooled objective still has enough freedom in directions the
    // regularizer is there to discourage. See heston_calibration.hpp.
    deriv::HestonCalibrationResult result = deriv::calibrate_heston(pooled, spot, r, guess, 0.001);

    std::printf("{\n");
    std::printf("  \"expiries_days\": [");
    for (std::size_t i = 0; i < expiry_report.size(); ++i) {
        std::printf("%s%.1f", i == 0 ? "" : ", ", expiry_report[i].first);
    }
    std::printf("],\n");
    std::printf("  \"contracts_per_expiry\": [");
    for (std::size_t i = 0; i < expiry_report.size(); ++i) {
        std::printf("%s%d", i == 0 ? "" : ", ", expiry_report[i].second);
    }
    std::printf("],\n");
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
