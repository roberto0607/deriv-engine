#include "deriv-engine/pde_solver.hpp"
#include <cmath>
#include <algorithm>

namespace deriv {

namespace {

double payoff(double S, double K, OptionType type) {
    if (type == OptionType::Call) {
        return std::max(S - K, 0.0);
    } else {
        return std::max(K - S, 0.0);
    }
}

}  // namespace

std::vector<double> solve_tridiagonal(const std::vector<double>& lower,
                                       const std::vector<double>& diag,
                                       const std::vector<double>& upper,
                                       const std::vector<double>& d) {
    std::size_t n = diag.size();
    std::vector<double> c_prime(n);  // modified super-diagonal
    std::vector<double> d_prime(n);  // modified right-hand side
    std::vector<double> x(n);

    // Forward sweep: eliminate the sub-diagonal, turning the system
    // into pure back-substitution — this is the Thomas algorithm,
    // the standard efficient method for tridiagonal systems (much
    // cheaper than general Gaussian elimination, which would waste
    // time on all the entries we already know are zero).
    c_prime[0] = upper[0] / diag[0];
    d_prime[0] = d[0] / diag[0];

    for (std::size_t i = 1; i < n; ++i) {
        double denom = diag[i] - lower[i] * c_prime[i - 1];
        c_prime[i] = upper[i] / denom;
        d_prime[i] = (d[i] - lower[i] * d_prime[i - 1]) / denom;
    }

    // Back substitution
    x[n - 1] = d_prime[n - 1];
    for (std::size_t i = n - 1; i-- > 0;) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    return x;
}

}  // namespace deriv