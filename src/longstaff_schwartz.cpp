#include "deriv-engine/longstaff_schwartz.hpp"
#include <cmath>
#include <random>
#include <vector>
#include <array>

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

// Fits Y ~ a0 + a1*X + a2*X^2 by least squares, using the "normal
// equations" approach: this reduces regression to solving a small
// 3x3 linear system, which we solve directly by Gaussian elimination.
// Returns {a0, a1, a2}.
std::array<double, 3> fit_quadratic(const std::vector<double>& X, const std::vector<double>& Y) {
    double n = static_cast<double>(X.size());
    double sx = 0, sx2 = 0, sx3 = 0, sx4 = 0;
    double sy = 0, sxy = 0, sx2y = 0;

    for (std::size_t i = 0; i < X.size(); ++i) {
        double x = X[i], y = Y[i];
        double x2 = x * x;
        sx += x;
        sx2 += x2;
        sx3 += x2 * x;
        sx4 += x2 * x2;
        sy += y;
        sxy += x * y;
        sx2y += x2 * y;
    }

    // Augmented matrix for [n, sx, sx2 | sy; sx, sx2, sx3 | sxy; sx2, sx3, sx4 | sx2y]
    double A[3][4] = {
        {n, sx, sx2, sy},
        {sx, sx2, sx3, sxy},
        {sx2, sx3, sx4, sx2y}
    };

    // Gaussian elimination with partial pivoting (swap rows to avoid
    // dividing by a near-zero pivot, which would blow up numerically).
    for (int col = 0; col < 3; ++col) {
        int pivot_row = col;
        for (int row = col + 1; row < 3; ++row) {
            if (std::abs(A[row][col]) > std::abs(A[pivot_row][col])) pivot_row = row;
        }
        std::swap(A[col], A[pivot_row]);

        for (int row = col + 1; row < 3; ++row) {
            double factor = A[row][col] / A[col][col];
            for (int k = col; k < 4; ++k) {
                A[row][k] -= factor * A[col][k];
            }
        }
    }

    // Back-substitution
    std::array<double, 3> coeffs;
    for (int row = 2; row >= 0; --row) {
        double sum = A[row][3];
        for (int col = row + 1; col < 3; ++col) {
            sum -= A[row][col] * coeffs[col];
        }
        coeffs[row] = sum / A[row][row];
    }

    return coeffs;
}

}  // namespace deriv