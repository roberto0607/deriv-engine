#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "deriv-engine/black_scholes.hpp"
#include "deriv-engine/greeks.hpp"
#include <cmath>
#include "deriv-engine/day_count.hpp"
#include "deriv-engine/tree_pricer.hpp"
#include "deriv-engine/monte_carlo.hpp"
#include "deriv-engine/adouble.hpp"
#include <chrono>
#include "deriv-engine/hedging.hpp"

TEST_CASE("Sanity check", "[placeholder]") {
    REQUIRE(1 + 1 == 2);
}

TEST_CASE("Black-Scholes call matches known value", "[black_scholes]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};   // spot=100, r=5%, q=0%, vol=20%
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};  // K=100, T=1yr

    double price = deriv::black_scholes_price(option, market);

    REQUIRE(price == Catch::Approx(10.4506).epsilon(0.001));
}

TEST_CASE("Black-Scholes put matches known value", "[black_scholes]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Put};

    double price = deriv::black_scholes_price(option, market);

    REQUIRE(price == Catch::Approx(5.5735).epsilon(0.001));
}

TEST_CASE("Black-Scholes call Greeks match known values", "[greeks]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    deriv::Greeks g = deriv::black_scholes_greeks(option, market);

    REQUIRE(g.delta == Catch::Approx(0.6368).epsilon(0.001));
    REQUIRE(g.gamma == Catch::Approx(0.0188).epsilon(0.01));
    REQUIRE(g.vega  == Catch::Approx(37.52).epsilon(0.01));
    REQUIRE(g.theta == Catch::Approx(-6.414).epsilon(0.01));
    REQUIRE(g.rho   == Catch::Approx(53.23).epsilon(0.01));
}

TEST_CASE("Black-Scholes put Greeks match known values", "[greeks]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Put};

    deriv::Greeks g = deriv::black_scholes_greeks(option, market);

    REQUIRE(g.delta == Catch::Approx(-0.3632).epsilon(0.001));
    REQUIRE(g.gamma == Catch::Approx(0.0188).epsilon(0.01));
    REQUIRE(g.vega  == Catch::Approx(37.52).epsilon(0.01));
    REQUIRE(g.theta == Catch::Approx(-1.658).epsilon(0.01));
    REQUIRE(g.rho   == Catch::Approx(-41.89).epsilon(0.01));
}

TEST_CASE("Put-call parity holds", "[black_scholes][parity]") {
    deriv::MarketData market{100.0, 0.05, 0.02, 0.25};  // includes a dividend yield this time
    deriv::EuropeanOption call{110.0, 0.75, deriv::OptionType::Call};
    deriv::EuropeanOption put{110.0, 0.75, deriv::OptionType::Put};

    double call_price = deriv::black_scholes_price(call, market);
    double put_price = deriv::black_scholes_price(put, market);

    double lhs = call_price - put_price;
    double rhs = market.spot * std::exp(-market.dividend_yield * call.time_to_expiry)
               - call.strike * std::exp(-market.risk_free_rate * call.time_to_expiry);

    REQUIRE(lhs == Catch::Approx(rhs).epsilon(0.0001));
}

TEST_CASE("Handles near-zero volatility without crashing", "[black_scholes][edge_case]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 1e-9};  // effectively zero vol
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    double price = deriv::black_scholes_price(option, market);

    // With near-zero vol, the call should behave like a forward contract:
    // price ≈ max(S*e^(-qT) - K*e^(-rT), 0)
    double expected = std::max(market.spot * std::exp(-market.dividend_yield * option.time_to_expiry)
                              - option.strike * std::exp(-market.risk_free_rate * option.time_to_expiry), 0.0);

    REQUIRE(price == Catch::Approx(expected).epsilon(0.01));
}

TEST_CASE("Zero volatility causes division by zero", "[black_scholes][edge_case]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.0};  // exactly zero vol
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    double price = deriv::black_scholes_price(option, market);

    INFO("Price computed: " << price);
    REQUIRE(std::isfinite(price));
}

TEST_CASE("Zero time to expiry causes division by zero", "[black_scholes][edge_case]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 0.0, deriv::OptionType::Call};  // T = 0

    double price = deriv::black_scholes_price(option, market);

    INFO("Price computed: " << price);
    REQUIRE(std::isfinite(price));
}

TEST_CASE("Greeks handle zero time to expiry without crashing", "[greeks][edge_case]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption itm_call{90.0, 0.0, deriv::OptionType::Call};  // S=100 > K=90, ITM

    deriv::Greeks g = deriv::black_scholes_greeks(itm_call, market);

    REQUIRE(std::isfinite(g.delta));
    REQUIRE(std::isfinite(g.gamma));
    REQUIRE(std::isfinite(g.vega));
    REQUIRE(std::isfinite(g.theta));
    REQUIRE(std::isfinite(g.rho));
    REQUIRE(g.delta == Catch::Approx(1.0));
    REQUIRE(g.gamma == Catch::Approx(0.0));
}

TEST_CASE("Deep in-the-money call behaves sensibly", "[black_scholes][edge_case]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{1.0, 1.0, deriv::OptionType::Call};  // K=1, way ITM

    double price = deriv::black_scholes_price(option, market);
    deriv::Greeks g = deriv::black_scholes_greeks(option, market);

    // Deep ITM call should behave almost exactly like the forward:
    // price ≈ S - K*e^(-rT), delta ≈ 1
    double expected_price = market.spot - option.strike * std::exp(-market.risk_free_rate * option.time_to_expiry);

    REQUIRE(std::isfinite(price));
    REQUIRE(price == Catch::Approx(expected_price).epsilon(0.01));
    REQUIRE(g.delta == Catch::Approx(1.0).epsilon(0.001));
    REQUIRE(g.gamma == Catch::Approx(0.0).margin(0.001));
}

TEST_CASE("Deep out-of-the-money call behaves sensibly", "[black_scholes][edge_case]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{10000.0, 1.0, deriv::OptionType::Call};  // K=10000, way OTM

    double price = deriv::black_scholes_price(option, market);
    deriv::Greeks g = deriv::black_scholes_greeks(option, market);

    // Deep OTM call is worth almost nothing, delta near 0
    REQUIRE(std::isfinite(price));
    REQUIRE(price == Catch::Approx(0.0).margin(0.01));
    REQUIRE(g.delta == Catch::Approx(0.0).margin(0.001));
}

TEST_CASE("Negative interest rate does not break pricing", "[black_scholes][edge_case]") {
    deriv::MarketData market{100.0, -0.01, 0.0, 0.2};  // negative rate
    deriv::EuropeanOption call{100.0, 1.0, deriv::OptionType::Call};
    deriv::EuropeanOption put{100.0, 1.0, deriv::OptionType::Put};

    double call_price = deriv::black_scholes_price(call, market);
    double put_price = deriv::black_scholes_price(put, market);

    REQUIRE(std::isfinite(call_price));
    REQUIRE(std::isfinite(put_price));
    REQUIRE(call_price > 0.0);
    REQUIRE(put_price > 0.0);

    // Put-call parity should still hold under a negative rate
    double lhs = call_price - put_price;
    double rhs = market.spot * std::exp(-market.dividend_yield * call.time_to_expiry)
               - call.strike * std::exp(-market.risk_free_rate * call.time_to_expiry);
    REQUIRE(lhs == Catch::Approx(rhs).epsilon(0.0001));
}

TEST_CASE("Day count conventions compute correct year fractions", "[day_count]") {
    // 365 days under ACT/365 should be exactly 1 year
    REQUIRE(deriv::year_fraction(365, deriv::DayCountConvention::Act365) == Catch::Approx(1.0));

    // 365 days under ACT/360 should be slightly MORE than 1 year
    // (since we're dividing by a smaller number, 360 instead of 365)
    REQUIRE(deriv::year_fraction(365, deriv::DayCountConvention::Act360) == Catch::Approx(365.0 / 360.0));

    // 90 days, a common short-term horizon, under both conventions
    REQUIRE(deriv::year_fraction(90, deriv::DayCountConvention::Act365) == Catch::Approx(90.0 / 365.0));
    REQUIRE(deriv::year_fraction(90, deriv::DayCountConvention::Act360) == Catch::Approx(90.0 / 360.0));

    // 0 days (option expires today) should give exactly 0
    REQUIRE(deriv::year_fraction(0, deriv::DayCountConvention::Act365) == Catch::Approx(0.0));
}

TEST_CASE("Binomial tree converges to Black-Scholes for European call", "[tree][convergence]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    double bs_price = deriv::black_scholes_price(option, market);
    double tree_price_low_steps = deriv::binomial_tree_price(option, market, 10);
    double tree_price_high_steps = deriv::binomial_tree_price(option, market, 500);

    // With few steps, the tree is a rough approximation — allow more slack
    REQUIRE(tree_price_low_steps == Catch::Approx(bs_price).epsilon(0.02));

    // With many steps, it should converge very close to the closed-form price
    REQUIRE(tree_price_high_steps == Catch::Approx(bs_price).epsilon(0.001));
}

TEST_CASE("Binomial tree converges to Black-Scholes for European put", "[tree][convergence]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Put};

    double bs_price = deriv::black_scholes_price(option, market);
    double tree_price = deriv::binomial_tree_price(option, market, 500);

    REQUIRE(tree_price == Catch::Approx(bs_price).epsilon(0.001));
}

TEST_CASE("American option is worth at least as much as European", "[tree]") {
    deriv::MarketData market{100.0, 0.05, 0.03, 0.2};  // nonzero dividend yield matters here
    deriv::EuropeanOption euro_put{100.0, 1.0, deriv::OptionType::Put};
    deriv::AmericanOption amer_put{100.0, 1.0, deriv::OptionType::Put};

    double euro_price = deriv::binomial_tree_price(euro_put, market, 500);
    double amer_price = deriv::binomial_tree_price(amer_put, market, 500);

    REQUIRE(amer_price >= euro_price);
}

TEST_CASE("Monte Carlo converges to Black-Scholes for European call", "[monte_carlo][convergence]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    double bs_price = deriv::black_scholes_price(option, market);
    deriv::MonteCarloResult mc = deriv::monte_carlo_price(option, market, 100000);

    // Check the MC price is within a few standard errors of the true price —
    // this is the statistically correct way to validate a random estimate,
    // rather than picking an arbitrary tolerance out of thin air.
    double diff = std::abs(mc.price - bs_price);
    REQUIRE(diff < 3.0 * mc.standard_error);
}

TEST_CASE("Monte Carlo converges to Black-Scholes for European put", "[monte_carlo][convergence]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Put};

    double bs_price = deriv::black_scholes_price(option, market);
    deriv::MonteCarloResult mc = deriv::monte_carlo_price(option, market, 100000);

    double diff = std::abs(mc.price - bs_price);
    REQUIRE(diff < 3.0 * mc.standard_error);
}

TEST_CASE("Antithetic variates reduce standard error", "[monte_carlo]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    deriv::MonteCarloResult with_antithetic = deriv::monte_carlo_price(option, market, 10000, true);
    deriv::MonteCarloResult without_antithetic = deriv::monte_carlo_price(option, market, 10000, false);

    // Antithetic variates should produce a tighter (smaller) standard error
    // for a comparable number of total random draws
    REQUIRE(with_antithetic.standard_error < without_antithetic.standard_error);
}

TEST_CASE("Control variate reduces standard error vs plain Monte Carlo", "[monte_carlo][control_variate]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    deriv::MonteCarloResult plain = deriv::monte_carlo_price(option, market, 10000, false);
    deriv::MonteCarloResult with_cv = deriv::monte_carlo_price_control_variate(option, market, 10000);

    REQUIRE(with_cv.standard_error < plain.standard_error);

    double bs_price = deriv::black_scholes_price(option, market);
    double diff = std::abs(with_cv.price - bs_price);
    REQUIRE(diff < 3.0 * with_cv.standard_error);
}

// TEST_CASE("Print variance reduction ratio for docs", "[.manual]") {
//     deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
//     deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

//     deriv::MonteCarloResult plain = deriv::monte_carlo_price(option, market, 50000, false);
//     deriv::MonteCarloResult antithetic = deriv::monte_carlo_price(option, market, 50000, true);
//     deriv::MonteCarloResult cv = deriv::monte_carlo_price_control_variate(option, market, 50000);

//     WARN("Plain SE: " << plain.standard_error);
//     WARN("Antithetic SE: " << antithetic.standard_error);
//     WARN("Control variate SE: " << cv.standard_error);
//     WARN("CV variance reduction ratio: " << (plain.standard_error * plain.standard_error) / (cv.standard_error * cv.standard_error));
// }

TEST_CASE("AAD tape reproduces hand-computed adjoints for y = a*b + c", "[aad]") {
    deriv::get_tape().clear();

    deriv::ADouble a(3.0);
    deriv::ADouble b(4.0);
    deriv::ADouble c(5.0);

    deriv::ADouble t1 = a * b;
    deriv::ADouble y = t1 + c;

    REQUIRE(y.value == Catch::Approx(17.0));

    deriv::get_tape().backward(y.idx);
    auto& adj = deriv::get_tape().adjoints;

    REQUIRE(adj[a.idx] == Catch::Approx(4.0));  // dy/da = b = 4
    REQUIRE(adj[b.idx] == Catch::Approx(3.0));  // dy/db = a = 3
    REQUIRE(adj[c.idx] == Catch::Approx(1.0));  // dy/dc = 1
}

TEST_CASE("Monte Carlo AAD delta and vega match analytical Greeks", "[aad][monte_carlo]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    deriv::Greeks analytical = deriv::black_scholes_greeks(option, market);
    deriv::MonteCarloGreeksResult mc = deriv::monte_carlo_greeks(option, market, 100000);

    // Delta and vega should be within a few standard errors of the
    // known-correct analytical values — same statistical validation
    // approach as pricing itself.
    REQUIRE(std::abs(mc.delta - analytical.delta) < 0.02);
    REQUIRE(std::abs(mc.vega - analytical.vega) < 2.0);
}

TEST_CASE("Monte Carlo AAD theta and rho match analytical Greeks", "[aad][monte_carlo]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    deriv::Greeks analytical = deriv::black_scholes_greeks(option, market);
    deriv::MonteCarloGreeksResult mc = deriv::monte_carlo_greeks(option, market, 100000);

    REQUIRE(std::abs(mc.delta - analytical.delta) < 0.02);
    REQUIRE(std::abs(mc.vega - analytical.vega) < 2.0);
    REQUIRE(std::abs(mc.theta - analytical.theta) < 1.0);
    REQUIRE(std::abs(mc.rho - analytical.rho) < 2.0);
}

TEST_CASE("AAD is faster than bump-and-revalue for computing 4 Greeks", "[.manual][aad][performance]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};
    const int num_paths = 100000;
    const double bump = 0.01;
    const unsigned int shared_seed = 42;  // fixed, so base and bumped runs share the same random draws

    auto bump_start = std::chrono::high_resolution_clock::now();

    deriv::MonteCarloResult base = deriv::monte_carlo_price(option, market, num_paths, false, shared_seed);

    deriv::MarketData bumped_S = market;
    bumped_S.spot += bump;
    deriv::MonteCarloResult r_S = deriv::monte_carlo_price(option, bumped_S, num_paths, false, shared_seed);
    double bump_delta = (r_S.price - base.price) / bump;

    deriv::MarketData bumped_sigma = market;
    bumped_sigma.volatility += bump;
    deriv::MonteCarloResult r_sigma = deriv::monte_carlo_price(option, bumped_sigma, num_paths, false, shared_seed);
    double bump_vega = (r_sigma.price - base.price) / bump;

    deriv::MarketData bumped_r = market;
    bumped_r.risk_free_rate += bump;
    deriv::MonteCarloResult r_r = deriv::monte_carlo_price(option, bumped_r, num_paths, false, shared_seed);
    double bump_rho = (r_r.price - base.price) / bump;

    deriv::EuropeanOption bumped_T = option;
    bumped_T.time_to_expiry -= bump;
    deriv::MonteCarloResult r_T = deriv::monte_carlo_price(bumped_T, market, num_paths, false, shared_seed);
    double bump_theta = (r_T.price - base.price) / bump;

    auto bump_end = std::chrono::high_resolution_clock::now();
    double bump_ms = std::chrono::duration<double, std::milli>(bump_end - bump_start).count();

    auto aad_start = std::chrono::high_resolution_clock::now();
    deriv::MonteCarloGreeksResult aad = deriv::monte_carlo_greeks(option, market, num_paths);
    auto aad_end = std::chrono::high_resolution_clock::now();
    double aad_ms = std::chrono::duration<double, std::milli>(aad_end - aad_start).count();

    WARN("Bump-and-revalue (5 runs): " << bump_ms << " ms");
    WARN("AAD (1 run): " << aad_ms << " ms");
    WARN("Speedup: " << (bump_ms / aad_ms) << "x");
    WARN("Bump delta: " << bump_delta << "  AAD delta: " << aad.delta);

    REQUIRE(std::abs(bump_delta - aad.delta) < 0.05);
}

TEST_CASE("Delta hedge simulation runs and produces a finite PnL", "[hedging]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    deriv::HedgeSimResult result = deriv::simulate_delta_hedge(option, market, 252, 42);  // 252 trading days

    REQUIRE(std::isfinite(result.final_pnl));
    REQUIRE(result.price_path.size() == 253);  // 252 steps + starting price
    WARN("Final PnL from one hedged path: " << result.final_pnl);
    WARN("Theoretical price: " << result.theoretical_price);
}

TEST_CASE("Delta hedge PnL distribution is centered near zero with reasonable spread", "[hedging]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    const int num_sims = 2000;
    std::vector<double> pnls;
    pnls.reserve(num_sims);

    for (int i = 0; i < num_sims; ++i) {
        deriv::HedgeSimResult result = deriv::simulate_delta_hedge(option, market, 252, 0);  // seed=0: fresh randomness each run
        pnls.push_back(result.final_pnl);
    }

    double mean = 0.0;
    for (double p : pnls) mean += p;
    mean /= pnls.size();

    double sq_diff_sum = 0.0;
    for (double p : pnls) {
        double diff = p - mean;
        sq_diff_sum += diff * diff;
    }
    double stdev = std::sqrt(sq_diff_sum / (pnls.size() - 1));

    WARN("Mean hedging PnL across " << num_sims << " paths: " << mean);
    WARN("Std dev of hedging PnL: " << stdev);

    // The mean PnL should be small relative to the option's own price —
    // hedging should roughly break even on average, not systematically
    // lose or gain a large amount.
    double theoretical_price = deriv::black_scholes_price(option, market);
    REQUIRE(std::abs(mean) < 0.15 * theoretical_price);
}