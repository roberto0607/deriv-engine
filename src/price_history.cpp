#include "deriv-engine/price_history.hpp"
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>

namespace deriv {

namespace {

// The actual "date,price" CSV parser -- shared by load_price_history()
// (reads a file) and load_price_history_from_string() (reads an
// in-memory string, for the WASM demo) so there's exactly one place
// this format is understood, not two copies that could drift apart.
std::vector<PricePoint> parse_price_history(std::istream& in) {
    std::vector<PricePoint> points;
    std::string line;

    // Header row, skipped unconditionally -- column names aren't
    // validated since only position ("date,price") matters here.
    std::getline(in, line);

    while (std::getline(in, line)) {
        if (line.empty()) continue;

        size_t comma = line.find(',');
        if (comma == std::string::npos) {
            throw std::runtime_error("load_price_history: malformed row (no comma): " + line);
        }

        std::string date = line.substr(0, comma);
        std::string price_str = line.substr(comma + 1);

        double price;
        try {
            price = std::stod(price_str);
        } catch (const std::exception&) {
            throw std::runtime_error("load_price_history: malformed price field in row: " + line);
        }

        points.push_back(PricePoint{date, price});
    }

    return points;
}

}  // namespace

std::vector<PricePoint> load_price_history(const std::string& csv_path) {
    std::ifstream file(csv_path);
    if (!file.is_open()) {
        throw std::runtime_error("load_price_history: could not open " + csv_path);
    }
    return parse_price_history(file);
}

std::vector<PricePoint> load_price_history_from_string(const std::string& csv_text) {
    std::istringstream stream(csv_text);
    return parse_price_history(stream);
}

}  // namespace deriv
