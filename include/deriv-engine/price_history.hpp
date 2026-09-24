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

// Same parser, same validation, same "date,price" format -- just fed
// from an in-memory string instead of a file path. Exists for the WASM
// demo (wasm/wasm_bridge.cpp): the browser fetches
// docs/btc_price_history.csv itself (there's no filesystem for
// load_price_history() to open inside a WASM module), and hands the raw
// text straight to this function, so parsing is still done by this
// project's real C++ code -- not reimplemented in JavaScript.
std::vector<PricePoint> load_price_history_from_string(const std::string& csv_text);

}  // namespace deriv
