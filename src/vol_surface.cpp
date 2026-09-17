#include "deriv-engine/vol_surface.hpp"
#include "deriv-engine/implied_vol.hpp"
#include "deriv-engine/market_data.hpp"
#include <nlohmann/json.hpp>
#include <fstream>
#include <ctime>
#include <cctype>

namespace deriv {

namespace {

int month_from_abbreviation(const std::string& abbr) {
    static const std::vector<std::string> months = {
        "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
        "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
    };
    for (std::size_t i = 0; i < months.size(); ++i) {
        if (abbr == months[i]) return static_cast<int>(i) + 1;
    }
    return -1;  // unrecognized — caller should treat as a parse failure
}

// Parses Deribit's instrument naming convention, e.g.
// "BTC-30OCT26-90000-P" -> day=30, month=OCT, year=2026, strike=90000,
// type=Put. Day length varies (1 or 2 digits), which is why we scan
// leading digits rather than assuming a fixed-width format.
struct ParsedInstrument {
    bool ok = false;
    int day = 0, month = 0, year = 0;
    double strike = 0.0;
    OptionType type = OptionType::Call;
};

ParsedInstrument parse_instrument_name(const std::string& name) {
    ParsedInstrument result;

    std::vector<std::string> parts;
    std::string current;
    for (char c : name) {
        if (c == '-') {
            parts.push_back(current);
            current.clear();
        } else {
            current += c;
        }
    }
    parts.push_back(current);

    if (parts.size() != 4) return result;  // not a standard option instrument name

    const std::string& date_part = parts[1];
    std::size_t i = 0;
    std::string day_str;
    while (i < date_part.size() && std::isdigit(static_cast<unsigned char>(date_part[i]))) {
        day_str += date_part[i];
        ++i;
    }
    if (day_str.empty() || i + 3 > date_part.size()) return result;

    std::string month_str = date_part.substr(i, 3);
    std::string year_str = date_part.substr(i + 3);
    if (year_str.size() != 2) return result;

    result.day = std::stoi(day_str);
    result.month = month_from_abbreviation(month_str);
    result.year = 2000 + std::stoi(year_str);
    if (result.month == -1) return result;

    try {
        result.strike = std::stod(parts[2]);
    } catch (...) {
        return result;
    }

    if (parts[3] == "P") {
        result.type = OptionType::Put;
    } else if (parts[3] == "C") {
        result.type = OptionType::Call;
    } else {
        return result;
    }

    result.ok = true;
    return result;
}

double years_between(int y1, int m1, int d1, int y2, int m2, int d2) {
    std::tm tm1{}, tm2{};
    tm1.tm_year = y1 - 1900; tm1.tm_mon = m1 - 1; tm1.tm_mday = d1;
    tm2.tm_year = y2 - 1900; tm2.tm_mon = m2 - 1; tm2.tm_mday = d2;
    std::time_t t1 = std::mktime(&tm1);
    std::time_t t2 = std::mktime(&tm2);
    double days = std::difftime(t2, t1) / (60.0 * 60.0 * 24.0);
    return days / 365.0;  // ACT/365, consistent with Phase 1's day-count convention
}

}  // namespace

std::vector<VolSurfacePoint> build_vol_surface_from_snapshot(const std::string& json_path,
                                                               int today_year, int today_month, int today_day) {
    std::vector<VolSurfacePoint> points;

    std::ifstream file(json_path);
    nlohmann::json data;
    file >> data;

    for (const auto& entry : data["result"]) {
        std::string name = entry.value("instrument_name", "");
        ParsedInstrument parsed = parse_instrument_name(name);
        if (!parsed.ok) continue;

        double mark_price_btc = entry.value("mark_price", 0.0);
        double underlying_price = entry.value("underlying_price", 0.0);
        double mark_iv_percent = entry.value("mark_iv", 0.0);
        double r = entry.value("interest_rate", 0.0);

        if (mark_price_btc <= 0.0 || underlying_price <= 0.0) continue;

        double T = years_between(today_year, today_month, today_day, parsed.year, parsed.month, parsed.day);
        if (T <= 0.0) continue;  // already expired relative to our "today"

        double market_price_usd = mark_price_btc * underlying_price;

        EuropeanOption option{parsed.strike, T, parsed.type};
        MarketData market{underlying_price, r, 0.0, 0.0};  // vol placeholder, ignored by the solver

        ImpliedVolResult iv = implied_volatility(option, market, market_price_usd);

        VolSurfacePoint point;
        point.strike = parsed.strike;
        point.time_to_expiry = T;
        point.moneyness = parsed.strike / underlying_price;
        point.type = parsed.type;
        point.market_iv_reported = mark_iv_percent / 100.0;
        point.solved_iv = iv.vol;
        point.converged = iv.converged;
        points.push_back(point);
    }

    return points;
}

}  // namespace deriv