// Regenerates the real headline numbers shown on docs/index.html from the
// actual engine, run against the committed real Deribit snapshot -- the
// same computation as the "[.manual]" vol-surface test in tests/test_main.cpp.
// Never hand-type these numbers into the page; this program is the only
// source of truth for them. Emits JSON to stdout.

#include "deriv-engine/vol_surface.hpp"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

int main() {
    auto points = deriv::build_vol_surface_from_snapshot("data/btc_chain_snapshot.json", 2026, 9, 17);

    const int total = static_cast<int>(points.size());
    int converged_count = 0;
    double total_abs_diff = 0.0;
    std::vector<double> gaps;
    for (const auto& p : points) {
        if (p.converged) {
            converged_count++;
            double gap = std::abs(p.solved_iv - p.market_iv_reported);
            total_abs_diff += gap;
            gaps.push_back(gap);
        }
    }

    const double avg_diff = converged_count > 0 ? total_abs_diff / converged_count : 0.0;
    const double convergence_pct = total > 0 ? (100.0 * converged_count / total) : 0.0;

    // Median and 90th-percentile gap, alongside the mean. The mean alone
    // is skewed by a small tail of genuinely hard contracts (1-day-to-expiry,
    // deep ITM/OTM) where the exchange's own quoted IV is likely thin or
    // stale to begin with -- median is the more honest "typical contract"
    // number, and p90 makes the tail's size explicit rather than hiding it.
    double median_diff = 0.0;
    double p90_diff = 0.0;
    if (!gaps.empty()) {
        std::sort(gaps.begin(), gaps.end());
        median_diff = gaps[gaps.size() / 2];
        p90_diff = gaps[static_cast<std::size_t>(gaps.size() * 0.9)];
    }

    std::printf("{\n");
    std::printf("  \"total_contracts\": %d,\n", total);
    std::printf("  \"converged_contracts\": %d,\n", converged_count);
    std::printf("  \"convergence_pct\": %.1f,\n", convergence_pct);
    std::printf("  \"median_gap_pp\": %.2f,\n", median_diff * 100.0);
    std::printf("  \"avg_gap_pp\": %.1f,\n", avg_diff * 100.0);
    std::printf("  \"p90_gap_pp\": %.1f\n", p90_diff * 100.0);
    std::printf("}\n");

    return 0;
}
