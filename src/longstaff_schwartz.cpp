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

LSMResult longstaff_schwartz_price(const AmericanOption& option, const MarketData& market,
                                    int num_paths, int num_steps, unsigned int seed) {
    const double S0 = market.spot;
    const double K = option.strike;
    const double T = option.time_to_expiry;
    const double r = market.risk_free_rate;
    const double q = market.dividend_yield;
    const double sigma = market.volatility;
    const double dt = T / num_steps;
    const double discount_per_step = std::exp(-r * dt);

    std::mt19937 rng(seed == 0 ? std::random_device{}() : seed);
    std::normal_distribution<double> normal(0.0, 1.0);

    // Step 1: simulate every path forward, storing the FULL price
    // history — unlike European pricing, we need the price at every
    // intermediate date, not just at expiry.
    std::vector<std::vector<double>> paths(num_paths, std::vector<double>(num_steps + 1));
    for (int p = 0; p < num_paths; ++p) {
        paths[p][0] = S0;
        for (int step = 1; step <= num_steps; ++step) {
            double z = normal(rng);
            double drift = (r - q - 0.5 * sigma * sigma) * dt;
            double diffusion = sigma * std::sqrt(dt) * z;
            paths[p][step] = paths[p][step - 1] * std::exp(drift + diffusion);
        }
    }

    // cash_flow[p] tracks the payoff each path will actually realize,
    // exercise_step[p] tracks WHEN, so we discount correctly at the end.
    std::vector<double> cash_flow(num_paths);
    std::vector<int> exercise_step(num_paths, num_steps);
    for (int p = 0; p < num_paths; ++p) {
        cash_flow[p] = payoff(paths[p][num_steps], K, option.type);
    }

    // Step 2: walk BACKWARD from the second-to-last step to the first.
    for (int step = num_steps - 1; step >= 1; --step) {
        std::vector<double> X;
        std::vector<double> Y;
        std::vector<int> itm_path_indices;

        for (int p = 0; p < num_paths; ++p) {
            double S = paths[p][step];
            double exercise_value = payoff(S, K, option.type);

            // Only regress over IN-THE-MONEY paths — regressing over
            // paths where exercise isn't even a candidate wastes data
            // and distorts the fit where it matters most.
            if (exercise_value > 0.0) {
                int steps_ahead = exercise_step[p] - step;
                double discounted_future = cash_flow[p] * std::pow(discount_per_step, steps_ahead);
                X.push_back(S);
                Y.push_back(discounted_future);
                itm_path_indices.push_back(p);
            }
        }

        if (X.size() < 3) continue;  // not enough data to fit reliably

        auto coeffs = fit_quadratic(X, Y);

        for (std::size_t i = 0; i < itm_path_indices.size(); ++i) {
            int p = itm_path_indices[i];
            double S = X[i];
            double continuation_value = coeffs[0] + coeffs[1] * S + coeffs[2] * S * S;
            double exercise_value = payoff(S, K, option.type);

            if (exercise_value > continuation_value) {
                cash_flow[p] = exercise_value;
                exercise_step[p] = step;
            }
        }
    }

    // Step 3: discount every path's realized cash flow back to today
    // from whenever it was actually realized, then average.
    std::vector<double> discounted_payoffs(num_paths);
    for (int p = 0; p < num_paths; ++p) {
        discounted_payoffs[p] = cash_flow[p] * std::pow(discount_per_step, exercise_step[p]);
    }

    double sum = 0.0;
    for (double v : discounted_payoffs) sum += v;
    double price = sum / num_paths;

    double sq_diff_sum = 0.0;
    for (double v : discounted_payoffs) {
        double diff = v - price;
        sq_diff_sum += diff * diff;
    }
    double sample_variance = sq_diff_sum / (num_paths - 1);
    double standard_error = std::sqrt(sample_variance / num_paths);

    return LSMResult{price, standard_error};
}

}  // namespace deriv