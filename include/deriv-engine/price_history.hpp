#pragma once

#include <string>
#include <vector>

namespace deriv {

// One row of a daily price history: a date string ("YYYY-MM-DD", not
// parsed into a calendar type since nothing here needs calendar
// arithmetic -- only day-to-day spacing, which the caller gets for free
// from row order) and that day's closing price.
struct PricePoint {
    std::string date;
    double price;
};

// Loads a two-column "date,price" CSV (header row required, exact column
// names don't matter -- only column position does) into an ordered
// vector, oldest first. Throws std::runtime_error on a missing file or a
// row that doesn't parse as "text,number". Deliberately doesn't validate
// date continuity or ordering here -- see backtest.hpp's rolling-window
// engine, which is the actual consumer and is where gaps/duplicates would
// matter, for how (and whether) it depends on that.
std::vector<PricePoint> load_price_history(const std::string& csv_path);

}  // namespace deriv
