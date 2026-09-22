#pragma once

#include "instrument.hpp"
#include <string>
#include <vector>

namespace deriv {

struct VolSurfacePoint {
    double strike;
    double time_to_expiry;
    double moneyness;           // strike / underlying_price
    OptionType type;
    double market_iv_reported;  // Deribit's own mark_iv, as a fraction (not percent)
    double solved_iv;           // this engine's independently-solved implied vol
    bool converged;

    // Carried through for consumers (e.g. Heston calibration) that need
    // to reprice this exact contract rather than just compare vols.
    double underlying_price;
    double market_price_usd;
    double risk_free_rate;
};

// Parses a Deribit book-summary JSON snapshot and solves implied vol
// for every contract, independently of Deribit's own reported mark_iv.
// today_year/month/day pins the "as of" date used for time-to-expiry
// calculations, so results are reproducible from a frozen snapshot
// rather than depending on when the code happens to run.
std::vector<VolSurfacePoint> build_vol_surface_from_snapshot(const std::string& json_path,
                                                               int today_year, int today_month, int today_day);

}  // namespace deriv