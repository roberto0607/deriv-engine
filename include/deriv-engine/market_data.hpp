#pragma once

namespace deriv {

struct MarketData {
    double spot;
    double risk_free_rate;
    double dividend_yield;
    double volatility;
};

}  // namespace deriv