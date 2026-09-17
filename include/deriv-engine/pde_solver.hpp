#pragma once

#include "market_data.hpp"
#include "instrument.hpp"
#include <vector>

namespace deriv {

double crank_nicolson_price(const EuropeanOption& option, const MarketData& market,
                             int num_price_steps, int num_time_steps);

// A deliberately simple, fully EXPLICIT finite difference scheme —
// built only to demonstrate conditional stability, not for production
// pricing use. With too large a time step relative to the grid
// spacing, this scheme is known to diverge.
double explicit_fd_price(const EuropeanOption& option, const MarketData& market,
                          int num_price_steps, int num_time_steps);

// Exposed for testing. Solves a tridiagonal system Ax = d, where A has
// sub-diagonal `lower`, main diagonal `diag`, and super-diagonal
// `upper`. All four vectors must be the same length n; lower[0] and
// upper[n-1] are unused (there's no neighbor there) but kept for
// simple, uniform indexing.
std::vector<double> solve_tridiagonal(const std::vector<double>& lower,
                                       const std::vector<double>& diag,
                                       const std::vector<double>& upper,
                                       const std::vector<double>& d);

}  // namespace deriv