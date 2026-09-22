#pragma once

#include "market_data.hpp"
#include "instrument.hpp"

namespace deriv {

// Heston (1993) stochastic-volatility model:
//   dS = (r - q) S dt + sqrt(v) S dW1
//   dv = kappa (theta - v) dt + xi sqrt(v) dW2,   corr(dW1, dW2) = rho
struct HestonParams {
    double v0;     // initial variance
    double kappa;  // mean-reversion speed of variance
    double theta;  // long-run variance
    double xi;     // vol-of-vol
    double rho;    // spot-vol correlation
};

// Prices a European option under Heston via the COS method (Fang &
// Oosterlee, 2008), using the numerically stable "little trap"
// characteristic function (Albrecher, Mayer, Schoutens, Tistaert,
// 2007) to avoid branch-cut discontinuities in the principal complex
// log/sqrt that the original 1993 formula suffers from.
double heston_price(const EuropeanOption& option, const MarketData& market, const HestonParams& params);

}  // namespace deriv
