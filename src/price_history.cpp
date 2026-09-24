#include "deriv-engine/price_history.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace deriv {

std::vector<PricePoint> load_price_history(const std::string& csv_path) {
    std::ifstream file(csv_path);
    if (!file.is_open()) {
        throw std::runtime_error("load_price_history: could not open " + csv_path);
    }

    std::vector<PricePoint> points;
    std::string line;

    // Header row, skipped unconditionally -- column names aren't
    // validated since only position ("date,price") matters here.
    std::getline(file, line);

    while (std::getline(file, line)) {
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

}  // namespace deriv
