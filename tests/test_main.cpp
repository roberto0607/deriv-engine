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
#include "deriv-engine/adouble2.hpp"
#include "deriv-engine/tape2.hpp"
#include <chrono>
#include "deriv-engine/hedging.hpp"
#include "deriv-engine/longstaff_schwartz.hpp"
#include "deriv-engine/pde_solver.hpp"
#include "deriv-engine/implied_vol.hpp"
#include "deriv-engine/vol_surface.hpp"
#include "deriv-engine/heston.hpp"
#include "deriv-engine/heston_mc.hpp"
#include "deriv-engine/heston_calibration.hpp"
#include "deriv-engine/exotic_options.hpp"
#include "deriv-engine/bates.hpp"
#include "deriv-engine/bates_mc.hpp"
#include "deriv-engine/price_history.hpp"
#include "deriv-engine/backtest.hpp"
#include "deriv-engine/portfolio.hpp"
#include "deriv-engine/var.hpp"
#include <algorithm>
#include <vector>


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

TEST_CASE("Trinomial tree converges to Black-Scholes for European call", "[tree][trinomial][convergence]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    double bs_price = deriv::black_scholes_price(option, market);
    double tree_price_low_steps = deriv::trinomial_tree_price(option, market, 10);
    double tree_price_high_steps = deriv::trinomial_tree_price(option, market, 500);

    REQUIRE(tree_price_low_steps == Catch::Approx(bs_price).epsilon(0.02));
    REQUIRE(tree_price_high_steps == Catch::Approx(bs_price).epsilon(0.001));
}

TEST_CASE("Trinomial tree converges to Black-Scholes for European put", "[tree][trinomial][convergence]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Put};

    double bs_price = deriv::black_scholes_price(option, market);
    double tree_price = deriv::trinomial_tree_price(option, market, 500);

    REQUIRE(tree_price == Catch::Approx(bs_price).epsilon(0.001));
}

TEST_CASE("Trinomial American option is worth at least as much as European", "[tree][trinomial]") {
    deriv::MarketData market{100.0, 0.05, 0.03, 0.2};
    deriv::EuropeanOption euro_put{100.0, 1.0, deriv::OptionType::Put};
    deriv::AmericanOption amer_put{100.0, 1.0, deriv::OptionType::Put};

    double euro_price = deriv::trinomial_tree_price(euro_put, market, 500);
    double amer_price = deriv::trinomial_tree_price(amer_put, market, 500);

    REQUIRE(amer_price >= euro_price);
}

TEST_CASE("Binomial and trinomial trees agree with each other at high step counts",
          "[tree][trinomial][convergence]") {
    // Both trees are independent discretizations of the same continuous
    // process (2 branches vs. 3 per node, different u/d/p vs.
    // dx/pu/pm/pd derivations) -- at enough steps, both must converge
    // to the same place regardless of which discretization got there,
    // which is a stronger check than each one separately matching
    // Black-Scholes (that only rules out a shared error, and there
    // isn't one to share here: run_tree() and run_trinomial_tree() are
    // two independently-derived implementations).
    struct Scenario {
        double spot, strike, T, vol, r, q;
        deriv::OptionType type;
        const char* label;
    };
    std::vector<Scenario> scenarios = {
        {100.0, 100.0, 1.0, 0.2, 0.05, 0.0, deriv::OptionType::Call, "ATM call"},
        {100.0, 120.0, 0.5, 0.35, 0.03, 0.02, deriv::OptionType::Put, "OTM put, short-dated, with dividend"},
        {100000.0, 90000.0, 90.0 / 365.0, 0.8, 0.05, 0.0, deriv::OptionType::Call, "BTC-scale ITM call, high vol"},
    };

    for (const auto& s : scenarios) {
        INFO(s.label);
        deriv::MarketData market{s.spot, s.r, s.q, s.vol};
        deriv::EuropeanOption option{s.strike, s.T, s.type};

        double bino = deriv::binomial_tree_price(option, market, 800);
        double trino = deriv::trinomial_tree_price(option, market, 800);

        REQUIRE(trino == Catch::Approx(bino).epsilon(0.002));
    }
}

// Measured, not just claimed: does the trinomial tree actually converge
// faster than the binomial tree per step, as textbooks generally
// suggest (smoother/monotonic convergence from the extra branch)? Or is
// that a wash at this project's scale? See docs/phase-log.md Phase 10
// for the honest answer, the same way Phase 4 measured (rather than
// assumed) where AAD's speed advantage actually shows up.
TEST_CASE("Trinomial tree convergence rate vs. binomial tree, measured directly",
          "[tree][trinomial][.manual]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};
    double bs_price = deriv::black_scholes_price(option, market);

    for (int steps : {10, 25, 50, 100, 200, 400}) {
        double bino_err = std::abs(deriv::binomial_tree_price(option, market, steps) - bs_price);
        double trino_err = std::abs(deriv::trinomial_tree_price(option, market, steps) - bs_price);
        WARN("steps=" << steps << "  binomial error=" << bino_err << "  trinomial error=" << trino_err);
    }
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

// --- Second-order AAD (gamma) --------------------------------------
//
// Phase 4 deferred gamma, on the grounds that it needs "a tape-of-
// tapes structure to differentiate the backward pass itself" — a
// separate piece of machinery from the first-order Tape/ADouble above.
// Tape2/ADouble2/Dual (tape2.hpp/adouble2.hpp/dual.hpp) are that
// machinery: forward-over-reverse AD, verified below first against
// hand-computed toy examples (same discipline Phase 4 used for the
// original tape, "y = a*b + c"), then wired into two real uses with
// two very different outcomes — see docs/phase-log.md Phase 11.

TEST_CASE("Second-order AAD tape reproduces hand-computed derivatives", "[aad][gamma]") {
    // f(x) = x^3: f'=3x^2, f''=6x. At x=2: f=8, f'=12, f''=12.
    {
        deriv::get_tape2().clear();
        deriv::ADouble2 x(2.0, /*seed_dot=*/1.0);
        deriv::ADouble2 y = x * x * x;
        deriv::get_tape2().backward(y.idx);
        auto& adj = deriv::get_tape2().adjoints;

        REQUIRE(y.value.val == Catch::Approx(8.0));
        REQUIRE(adj[x.idx].val == Catch::Approx(12.0));  // f'
        REQUIRE(adj[x.idx].dot == Catch::Approx(12.0));  // f''
    }

    // f(x,y) = x^2*y + e^x, second derivative w.r.t. x only (y not
    // seeded) -- exercises a multivariate case, not just a single input.
    // d/dx = 2xy + e^x; d2/dx2 = 2y + e^x. At x=1.5, y=3.
    {
        deriv::get_tape2().clear();
        deriv::ADouble2 x(1.5, 1.0);
        deriv::ADouble2 y(3.0, 0.0);
        deriv::ADouble2 f = x * x * y + deriv::ad2_exp(x);
        deriv::get_tape2().backward(f.idx);
        auto& adj = deriv::get_tape2().adjoints;

        double e15 = std::exp(1.5);
        REQUIRE(adj[x.idx].val == Catch::Approx(2 * 1.5 * 3.0 + e15));
        REQUIRE(adj[x.idx].dot == Catch::Approx(2 * 3.0 + e15));
    }
}

TEST_CASE("Gamma via second-order AAD matches analytical Black-Scholes gamma", "[aad][gamma]") {
    // Applied to the Black-Scholes closed-form formula (a smooth
    // function of spot -- no kink), forward-over-reverse AD should
    // recover gamma to near machine precision, not just within a
    // statistical tolerance. This is the WORKING half of the gamma-via-
    // AAD story; see the next test for the other half.
    struct Scenario {
        deriv::MarketData market;
        deriv::EuropeanOption option;
        const char* label;
    };
    std::vector<Scenario> scenarios = {
        {{100.0, 0.05, 0.0, 0.20}, {100.0, 1.0, deriv::OptionType::Call}, "ATM call"},
        {{100.0, 0.05, 0.0, 0.20}, {100.0, 1.0, deriv::OptionType::Put}, "ATM put"},
        {{100.0, 0.05, 0.02, 0.45}, {130.0, 0.5, deriv::OptionType::Call}, "OTM call, high vol, dividend"},
        {{100.0, 0.03, 0.0, 0.25}, {150.0, 2.0, deriv::OptionType::Put}, "Deep ITM put, long-dated"},
        {{100000.0, 0.05, 0.0, 0.80}, {105000.0, 90.0 / 365.0, deriv::OptionType::Call}, "BTC-scale short-dated"},
        {{100.0, 0.05, 0.0, 0.2}, {100.0, 1.0 / 365.0, deriv::OptionType::Call}, "Near expiry"},
    };

    for (const auto& s : scenarios) {
        INFO(s.label);
        deriv::Greeks analytical = deriv::black_scholes_greeks(s.option, s.market);
        deriv::GammaAD2Result ad2 = deriv::black_scholes_greeks_via_ad2(s.option, s.market);
        double plain_price = deriv::black_scholes_price(s.option, s.market);

        REQUIRE(ad2.price == Catch::Approx(plain_price).epsilon(1e-9));
        REQUIRE(ad2.delta == Catch::Approx(analytical.delta).epsilon(1e-9));
        REQUIRE(ad2.gamma == Catch::Approx(analytical.gamma).epsilon(1e-9));
    }
}

TEST_CASE("Monte Carlo pathwise AAD gamma is exactly zero for a kinked payoff",
          "[aad][gamma][monte_carlo]") {
    // The other half of the story, and the real reason Phase 4 called
    // gamma a "genuinely larger sub-project": applying the exact same
    // second-order tape to a Monte Carlo path's payoff instead of the
    // smooth Black-Scholes formula does NOT recover Black-Scholes'
    // analytical gamma. It gives exactly zero, every time, regardless
    // of path count -- this is not noise or a bug in Tape2 (verified
    // correct above), it's a real, well-known limitation of pathwise
    // differentiation: each simulated path's payoff (max(S_T-K,0), or
    // the put equivalent) is PIECEWISE LINEAR in spot -- either exactly
    // S_T-K or exactly 0 -- and the second derivative of a piecewise
    // linear function is identically zero everywhere except exactly at
    // the kink (probability zero under continuous GBM). A working Monte
    // Carlo gamma estimator needs a fundamentally different technique
    // (e.g. the likelihood ratio / score-function method, or a smoothed
    // payoff) -- deliberately out of scope here; this test documents
    // WHY, with a real number, rather than silently shipping a gamma
    // that's always zero without explaining it.
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    deriv::MonteCarloGammaResult mc = deriv::monte_carlo_gamma(option, market, 50000, /*seed=*/7);

    REQUIRE(mc.gamma == 0.0);

    // Delta, however, comes from the SAME kind of pathwise
    // differentiation and works fine (the payoff is continuous, even
    // though its second derivative isn't) -- cross-checked against
    // monte_carlo_greeks' independently-computed (first-order tape)
    // delta, and against the analytical value, confirming this
    // function's own machinery is otherwise correct and the zero
    // gamma above isn't a sign that something else is broken.
    deriv::Greeks analytical = deriv::black_scholes_greeks(option, market);
    REQUIRE(std::abs(mc.delta - analytical.delta) < 0.02);
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

TEST_CASE("More frequent rebalancing reduces hedging error", "[.manual][hedging]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    auto run_experiment = [&](int num_rebalances) {
        const int num_sims = 2000;
        std::vector<double> pnls;
        pnls.reserve(num_sims);
        for (int i = 0; i < num_sims; ++i) {
            deriv::HedgeSimResult result = deriv::simulate_delta_hedge(option, market, num_rebalances, 0);
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
        return std::sqrt(sq_diff_sum / (pnls.size() - 1));
    };

    double stdev_monthly = run_experiment(12);
    double stdev_weekly = run_experiment(52);
    double stdev_daily = run_experiment(252);

    WARN("Hedging error stdev — monthly (12 rebalances): " << stdev_monthly);
    WARN("Hedging error stdev — weekly (52 rebalances): " << stdev_weekly);
    WARN("Hedging error stdev — daily (252 rebalances): " << stdev_daily);

    // More frequent rebalancing should shrink the hedging error —
    // this is the direct, quantified link back to Black-Scholes's
    // continuous-hedging assumption from Phase 1.
    REQUIRE(stdev_weekly < stdev_monthly);
    REQUIRE(stdev_daily < stdev_weekly);
}

TEST_CASE("Quadratic regression exactly recovers known coefficients with zero noise", "[lsm][regression]") {
    // Generate points that lie EXACTLY on y = 2 + 3x + 0.5x^2
    std::vector<double> X = {-2.0, -1.0, 0.0, 1.0, 2.0, 3.0};
    std::vector<double> Y;
    for (double x : X) {
        Y.push_back(2.0 + 3.0 * x + 0.5 * x * x);
    }

    auto coeffs = deriv::fit_quadratic(X, Y);

    REQUIRE(coeffs[0] == Catch::Approx(2.0).epsilon(0.0001));
    REQUIRE(coeffs[1] == Catch::Approx(3.0).epsilon(0.0001));
    REQUIRE(coeffs[2] == Catch::Approx(0.5).epsilon(0.0001));
}

TEST_CASE("Longstaff-Schwartz converges to tree price for American put", "[lsm][convergence]") {
    deriv::MarketData market{100.0, 0.05, 0.03, 0.2};  // nonzero dividend, matters for early exercise
    deriv::AmericanOption option{100.0, 1.0, deriv::OptionType::Put};

    double tree_price = deriv::binomial_tree_price(option, market, 500);
    deriv::LSMResult lsm = deriv::longstaff_schwartz_price(option, market, 50000, 50, 42);

    double diff = std::abs(lsm.price - tree_price);
    REQUIRE(diff < 3.0 * lsm.standard_error);
}

TEST_CASE("Longstaff-Schwartz converges to tree price for American call with dividends", "[lsm][convergence]") {
    // A call with a meaningful dividend yield — early exercise genuinely
    // matters here (unlike a non-dividend call, where it never would),
    // so this exercises the LSM decision logic in a different, equally
    // real scenario from the put case already tested.
    deriv::MarketData market{100.0, 0.03, 0.05, 0.25};
    deriv::AmericanOption option{95.0, 0.5, deriv::OptionType::Call};

    double tree_price = deriv::binomial_tree_price(option, market, 500);
    deriv::LSMResult lsm = deriv::longstaff_schwartz_price(option, market, 50000, 50, 123);

    double diff = std::abs(lsm.price - tree_price);
    REQUIRE(diff < 3.0 * lsm.standard_error);
}

TEST_CASE("Tridiagonal solver matches a hand-verified known solution", "[pde]") {
    // System:  2x0 -  x1        = 1
    //          -x0 + 2x1 -  x2  = 0
    //               -x1 + 2x2  = 1
    // Hand-solvable: this is a standard small tridiagonal test case,
    // with known exact solution x = [1, 1, 1].
    std::vector<double> lower = {0, -1, -1};   // lower[0] unused
    std::vector<double> diag  = {2, 2, 2};
    std::vector<double> upper = {-1, -1, 0};   // upper[2] unused
    std::vector<double> d     = {1, 0, 1};

    auto x = deriv::solve_tridiagonal(lower, diag, upper, d);

    REQUIRE(x[0] == Catch::Approx(1.0));
    REQUIRE(x[1] == Catch::Approx(1.0));
    REQUIRE(x[2] == Catch::Approx(1.0));
}

TEST_CASE("Crank-Nicolson converges to Black-Scholes for European call", "[pde][convergence]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    double bs_price = deriv::black_scholes_price(option, market);
    double pde_price = deriv::crank_nicolson_price(option, market, 200, 200);

    REQUIRE(pde_price == Catch::Approx(bs_price).epsilon(0.01));
}

TEST_CASE("Crank-Nicolson converges to Black-Scholes for European put", "[pde][convergence]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Put};

    double bs_price = deriv::black_scholes_price(option, market);
    double pde_price = deriv::crank_nicolson_price(option, market, 200, 200);

    REQUIRE(pde_price == Catch::Approx(bs_price).epsilon(0.01));
}

TEST_CASE("Explicit scheme diverges with too-large time step, Crank-Nicolson stays stable", "[.manual][pde][stability]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};
    double bs_price = deriv::black_scholes_price(option, market);

    // Deliberately choose a FINE price grid but a COARSE time grid —
    // this specific combination is what breaks the stability condition
    // for the explicit scheme (its stability depends on dt being small
    // relative to dS^2, not just small in absolute terms).
    int fine_price_steps = 200;
    int coarse_time_steps = 20;

    double explicit_price = deriv::explicit_fd_price(option, market, fine_price_steps, coarse_time_steps);
    double cn_price = deriv::crank_nicolson_price(option, market, fine_price_steps, coarse_time_steps);

    WARN("True Black-Scholes price: " << bs_price);
    WARN("Explicit scheme price (unstable regime): " << explicit_price);
    WARN("Crank-Nicolson price (same grid): " << cn_price);

    // The explicit scheme should be wildly wrong (or even nan/inf) in
    // this regime; Crank-Nicolson should still be close to correct.
    REQUIRE(std::abs(cn_price - bs_price) < 1.0);
}

TEST_CASE("Implied vol solver recovers known volatility from its own price", "[implied_vol]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};  // true vol = 20%
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    double known_price = deriv::black_scholes_price(option, market);

    deriv::ImpliedVolResult result = deriv::implied_volatility(option, market, known_price);

    REQUIRE(result.converged);
    REQUIRE(result.vol == Catch::Approx(0.2).epsilon(0.0001));
}

TEST_CASE("Implied vol solver recovers known vol across moneyness, maturity, and vol regimes",
          "[implied_vol][property]") {
    // Each scenario: generate a price from a KNOWN vol via Black-Scholes
    // (ground truth, not market data), then check the solver recovers
    // that exact vol. Unlike the single ATM case below, this sweeps
    // moneyness, time-to-expiry, and vol level together, so a change
    // that only breaks convergence in one corner of the parameter space
    // (e.g. long-dated, or very high vol) can't slip through on the
    // strength of the ATM case alone. All scenarios here have healthy
    // vega (checked independently) -- they exercise the normal
    // Newton-Raphson path, not the bisection fallback (see the
    // dedicated bisection test below for that).
    struct Scenario {
        double spot, strike, T, true_vol;
        deriv::OptionType type;
        const char* label;
    };

    std::vector<Scenario> scenarios = {
        {100.0, 100.0, 1.0, 0.20, deriv::OptionType::Call, "ATM, 1yr, 20% vol (baseline)"},
        {100.0, 100.0, 1.0 / 365.0, 0.20, deriv::OptionType::Call, "ATM, 1-day, 20% vol"},
        {100.0, 130.0, 1.0, 0.60, deriv::OptionType::Call, "30% OTM call, 1yr, 60% vol"},
        {100.0, 70.0, 1.0, 0.60, deriv::OptionType::Put, "30% OTM put, 1yr, 60% vol"},
        {100.0, 100.0, 3.0, 0.05, deriv::OptionType::Call, "ATM, 3yr, very low 5% vol"},
        {100.0, 100.0, 3.0, 1.50, deriv::OptionType::Call, "ATM, 3yr, extreme 150% vol"},
        {100000.0, 120000.0, 90.0 / 365.0, 0.80, deriv::OptionType::Call, "BTC 20% OTM call, 90d, 80% vol"},
        {100000.0, 80000.0, 90.0 / 365.0, 0.80, deriv::OptionType::Put, "BTC 20% OTM put, 90d, 80% vol"},
    };

    for (const auto& s : scenarios) {
        INFO(s.label);

        deriv::MarketData market_with_vol{s.spot, 0.05, 0.0, s.true_vol};
        deriv::EuropeanOption option{s.strike, s.T, s.type};
        double known_price = deriv::black_scholes_price(option, market_with_vol);

        deriv::MarketData market_without_vol{s.spot, 0.05, 0.0, 0.0};
        deriv::ImpliedVolResult result = deriv::implied_volatility(option, market_without_vol, known_price);

        REQUIRE(result.converged);
        REQUIRE(result.vol == Catch::Approx(s.true_vol).epsilon(0.001));
    }
}

TEST_CASE("Implied vol solver's bisection fallback converges (but honestly can't recover the exact "
          "vol) when vega underflows",
          "[implied_vol][edge_case][bisection]") {
    // A deep ITM BTC call one day from expiry: the option is nowhere
    // near worthless (price ~= intrinsic value, $60,005) but vega
    // underflows to effectively zero (checked below, ~1.6e-182) because
    // there's almost no time value left to be sensitive to vol. This is
    // exactly the case the bisection fallback in implied_volatility()
    // exists for -- Newton-Raphson can't take a step here since it
    // would be dividing by ~0.
    //
    // This test originally asserted the fallback recovers the *exact*
    // known vol -- that assertion was wrong, and failed when first
    // written, catching a real finding: when vega has underflowed this
    // far, the price is IDENTICAL to double-precision for any vol
    // across an enormous range (verified: 1% to 200% vol all produce
    // the exact same price here). No solver, however good, can recover
    // a meaningful vol from price alone in this regime -- the price
    // curve is flatter than the solver's own 1e-6 tolerance can
    // resolve. So the honest, correct claim for bisection here is
    // narrower: it finds *a* vol that reproduces the target price
    // within tolerance (converged == true, matching the old behavior's
    // outright non-convergence), not *the* vol that produced it. This
    // mirrors exactly why the live site's own calibration footnote
    // flags deep ITM/OTM, near-expiry contracts as unreliable for IV.
    double spot = 100000.0;
    double strike = 40000.0;
    double T = 1.0 / 365.0;
    double true_vol = 0.60;

    deriv::MarketData market_with_vol{spot, 0.05, 0.0, true_vol};
    deriv::EuropeanOption option{strike, T, deriv::OptionType::Call};

    // Confirm the premise: vega really is small enough to force the
    // fallback path, so this test keeps meaning what it says even if
    // implied_vol.cpp's internal threshold changes later.
    deriv::Greeks g = deriv::black_scholes_greeks(option, market_with_vol);
    INFO("vega at this point: " << g.vega);
    REQUIRE(std::abs(g.vega) < 1e-8);

    double known_price = deriv::black_scholes_price(option, market_with_vol);
    deriv::MarketData market_without_vol{spot, 0.05, 0.0, 0.0};
    deriv::ImpliedVolResult result = deriv::implied_volatility(option, market_without_vol, known_price);

    REQUIRE(result.converged);
    REQUIRE(std::isfinite(result.vol));

    // The recovered vol reproduces the target price -- that's the part
    // bisection can actually guarantee here.
    deriv::MarketData market_at_recovered{spot, 0.05, 0.0, result.vol};
    double price_at_recovered = deriv::black_scholes_price(option, market_at_recovered);
    REQUIRE(price_at_recovered == Catch::Approx(known_price).margin(1e-5));
}

TEST_CASE("Implied vol solver handles deep OTM gracefully (near-zero vega)", "[implied_vol][edge_case]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption option{10000.0, 1.0, deriv::OptionType::Call};  // deep OTM, from Phase 1's own test

    double tiny_price = 0.0001;  // an almost-worthless option price

    deriv::ImpliedVolResult result = deriv::implied_volatility(option, market, tiny_price);

    // We don't require convergence here — the point is it doesn't
    // crash, hang, or produce nan/inf, even in a genuinely hard case.
    REQUIRE(std::isfinite(result.vol));
}

TEST_CASE("Implied vol solver on real Deribit BTC put data", "[.manual][implied_vol][real_data]") {
    // Real snapshot pulled from Deribit's public API, BTC-30OCT26-90000-P
    double spot = 76871.89;
    double strike = 90000.0;
    double mark_price_btc = 0.176;
    double mark_price_usd = mark_price_btc * spot;
    double T = 43.0 / 365.0;  // ~43 days to expiry, ACT/365
    double r = 0.0;           // Deribit reports interest_rate: 0.0 for this instrument

    deriv::MarketData market{spot, r, 0.0, 0.0};  // vol placeholder, ignored by the solver
    deriv::EuropeanOption option{strike, T, deriv::OptionType::Put};

    deriv::ImpliedVolResult result = deriv::implied_volatility(option, market, mark_price_usd);

    WARN("Mark price (USD): " << mark_price_usd);
    WARN("Your solved implied vol: " << result.vol * 100.0 << "%");
    WARN("Deribit's reported mark_iv: 34.36%");
    WARN("Converged: " << result.converged << ", iterations: " << result.iterations);

    REQUIRE(result.converged);
}

TEST_CASE("Vol surface builds from real Deribit snapshot and mostly converges", "[.manual][vol_surface][real_data]") {
    auto points = deriv::build_vol_surface_from_snapshot("data/btc_chain_snapshot.json", 2026, 9, 17);

    REQUIRE(points.size() > 50);  // sanity check: real chain should have many contracts

    int converged_count = 0;
    double total_abs_diff = 0.0;
    for (const auto& p : points) {
        if (p.converged) {
            converged_count++;
            total_abs_diff += std::abs(p.solved_iv - p.market_iv_reported);
        }
    }

    double avg_diff = total_abs_diff / converged_count;
    WARN("Total contracts parsed: " << points.size());
    WARN("Converged: " << converged_count);
    WARN("Average |solved_iv - market_iv| across converged points: " << avg_diff);

    REQUIRE(converged_count > points.size() / 2);  // most should converge
}

TEST_CASE("Vol surface calibration numbers are pinned against a regression",
          "[.manual][vol_surface][real_data][regression]") {
    // The committed snapshot is fixed and the solver is deterministic
    // (no RNG involved, unlike Monte Carlo), so these numbers should be
    // exactly reproducible. This test exists so that a future change to
    // implied_volatility() or build_vol_surface_from_snapshot() that
    // silently regresses calibration quality fails CI, instead of only
    // being caught by eyeballing the live site. If a real, intentional
    // improvement moves these numbers, update the pinned values here —
    // don't just loosen the tolerance.
    auto points = deriv::build_vol_surface_from_snapshot("data/btc_chain_snapshot.json", 2026, 9, 17);

    REQUIRE(points.size() == 904);

    int converged_count = 0;
    std::vector<double> gaps;
    for (const auto& p : points) {
        if (p.converged) {
            converged_count++;
            gaps.push_back(std::abs(p.solved_iv - p.market_iv_reported));
        }
    }
    std::sort(gaps.begin(), gaps.end());
    double median_gap = gaps[gaps.size() / 2];

    REQUIRE(converged_count == 875);
    REQUIRE(median_gap == Catch::Approx(0.0026).margin(0.0005));
}

// --- Heston COS-method pricer ---------------------------------------------
//
// heston.cpp's own header comments document two independent correctness
// checks for the characteristic function: this deterministic-variance
// limit (below) and a Monte Carlo cross-check of the actual Heston SDEs
// (not yet implemented -- tracked as a follow-up). This limit test is the
// stronger of the two for catching CF bugs: with xi -> 0 the variance
// process has (almost) no randomness, v(t) stays at v0 = theta for all t,
// and the Heston model must degenerate EXACTLY to Black-Scholes with
// sigma = sqrt(v0). Any sign error in g(u) or the wrong branch of d(u)
// in heston_log_return_cf shows up here as a real pricing discrepancy,
// not just a numerical-stability artifact.
TEST_CASE("Heston collapses to Black-Scholes in the deterministic-variance limit",
          "[heston][property]") {
    struct Scenario {
        double spot, strike, T, vol;
        deriv::OptionType type;
        const char* label;
    };

    std::vector<Scenario> scenarios = {
        {100.0, 100.0, 1.0, 0.20, deriv::OptionType::Call, "ATM call, 1yr, 20% vol"},
        {100.0, 100.0, 1.0, 0.20, deriv::OptionType::Put, "ATM put, 1yr, 20% vol"},
        {100.0, 130.0, 1.0, 0.35, deriv::OptionType::Call, "30% OTM call, 1yr, 35% vol"},
        {100.0, 70.0, 1.0, 0.35, deriv::OptionType::Put, "30% OTM put, 1yr, 35% vol"},
        {100.0, 100.0, 1.0 / 365.0, 0.20, deriv::OptionType::Call, "ATM call, 1-day, 20% vol"},
        {100.0, 100.0, 3.0, 0.60, deriv::OptionType::Call, "ATM call, 3yr, 60% vol"},
        {100000.0, 120000.0, 90.0 / 365.0, 0.80, deriv::OptionType::Call, "BTC 20% OTM call, 90d, 80% vol"},
        {100000.0, 80000.0, 90.0 / 365.0, 0.80, deriv::OptionType::Put, "BTC 20% OTM put, 90d, 80% vol"},
        // Extreme moneyness (strike = 3x spot): found via manual testing
        // of the live demo (spot=100000, strike=300000, 90d, 60% vol) to
        // expose a real precision bug -- deriving a tiny call price via
        // put-call parity from a put ~200,000x larger loses so much
        // absolute precision in double that the "collapse to Black-
        // Scholes" property broke by ~82% here specifically, even though
        // every other scenario above (max ~1.3x moneyness) passed. Root
        // cause was catastrophic cancellation in heston_log_return_cf's
        // log((1-g*E)/(1-g)) term, NOT anything long double could fix --
        // an earlier long-double attempt "fixed" this on x86 but was a
        // silent no-op on Apple Silicon, where long double == double
        // (see docs/phase-log.md for the full story). The real,
        // portable fix uses a complex log1p (clog1p() in heston.cpp)
        // to compute that log accurately regardless of platform; this
        // case is the regression test for that fix.
        {100000.0, 300000.0, 90.0 / 365.0, 0.60, deriv::OptionType::Call, "BTC deep OTM call, strike=3x spot, 90d, 60% vol"},
    };

    // xi is small but not literally zero: heston_log_return_cf divides by
    // xi^2 in a couple of intermediate quantities (cancelled analytically,
    // per the comments in heston.cpp), and a genuinely zero xi would test
    // that cancellation rather than the deterministic-variance limit
    // itself. 1e-6 puts the variance process's own randomness many orders
    // of magnitude below the pricing tolerance below while still
    // exercising the real code path.
    const double xi = 1e-6;

    for (const auto& s : scenarios) {
        INFO(s.label);

        deriv::HestonParams params{s.vol * s.vol, 2.0, s.vol * s.vol, xi, -0.5};
        deriv::MarketData market{s.spot, 0.05, 0.0, 0.0};  // vol unused by Heston pricer
        deriv::EuropeanOption option{s.strike, s.T, s.type};

        double heston = deriv::heston_price(option, market, params);

        deriv::MarketData bs_market{s.spot, 0.05, 0.0, s.vol};
        double bs = deriv::black_scholes_price(option, bs_market);

        REQUIRE(heston == Catch::Approx(bs).epsilon(0.001));
    }
}

TEST_CASE("Heston call and put prices satisfy put-call parity", "[heston]") {
    // Checked directly against Heston's own put-call parity relation
    // (independent of Black-Scholes), across parameter sets that actually
    // exercise stochastic vol (nonzero xi, nonzero rho) -- unlike the
    // deterministic-variance limit test above, which deliberately keeps
    // xi negligible.
    struct Params {
        deriv::HestonParams heston;
        double spot, strike, T, r, q;
        const char* label;
    };

    std::vector<Params> cases = {
        {{0.04, 2.0, 0.04, 0.5, -0.7}, 100.0, 100.0, 1.0, 0.05, 0.0, "ATM, moderate vol-of-vol, negative rho"},
        {{0.64, 1.5, 0.64, 0.9, -0.6}, 100000.0, 110000.0, 90.0 / 365.0, 0.05, 0.0, "BTC-like, high vol-of-vol"},
        {{0.09, 3.0, 0.16, 0.3, 0.4}, 100.0, 90.0, 2.0, 0.03, 0.02, "positive rho, v0 != theta, with dividend"},
    };

    for (const auto& c : cases) {
        INFO(c.label);

        deriv::MarketData market{c.spot, c.r, c.q, 0.0};
        deriv::EuropeanOption call{c.strike, c.T, deriv::OptionType::Call};
        deriv::EuropeanOption put{c.strike, c.T, deriv::OptionType::Put};

        double call_price = deriv::heston_price(call, market, c.heston);
        double put_price = deriv::heston_price(put, market, c.heston);

        double disc_spot = c.spot * std::exp(-c.q * c.T);
        double disc_strike = c.strike * std::exp(-c.r * c.T);

        REQUIRE(call_price - put_price == Catch::Approx(disc_spot - disc_strike).margin(1e-6));
    }
}

TEST_CASE("Heston call price increases monotonically with initial variance", "[heston][property]") {
    // A basic model-sanity check independent of any closed-form
    // reference: holding everything else fixed, more initial variance
    // should never make a European call cheaper. Not a tight numerical
    // check, but exactly the kind of bug (e.g. a sign error reversing
    // the CF's dependence on v0) that could otherwise slip past the
    // deterministic-variance and put-call-parity tests above, since
    // both of those only probe specific slices of parameter space.
    deriv::MarketData market{100.0, 0.05, 0.0, 0.0};
    deriv::EuropeanOption option{100.0, 1.0, deriv::OptionType::Call};

    std::vector<double> v0_values = {0.01, 0.04, 0.09, 0.16, 0.25};
    double prev_price = -1.0;
    for (double v0 : v0_values) {
        deriv::HestonParams params{v0, 2.0, 0.04, 0.5, -0.5};
        double price = deriv::heston_price(option, market, params);
        REQUIRE(price > prev_price);
        prev_price = price;
    }
}

TEST_CASE("Heston price reduces to intrinsic value as expiry approaches", "[heston][edge_case]") {
    deriv::HestonParams params{0.04, 2.0, 0.04, 0.5, -0.5};

    deriv::MarketData market{100.0, 0.05, 0.0, 0.0};
    deriv::EuropeanOption itm_call{80.0, 0.0, deriv::OptionType::Call};
    deriv::EuropeanOption otm_call{120.0, 0.0, deriv::OptionType::Call};
    deriv::EuropeanOption itm_put{120.0, 0.0, deriv::OptionType::Put};

    REQUIRE(deriv::heston_price(itm_call, market, params) == Catch::Approx(20.0).margin(1e-9));
    REQUIRE(deriv::heston_price(otm_call, market, params) == Catch::Approx(0.0).margin(1e-9));
    REQUIRE(deriv::heston_price(itm_put, market, params) == Catch::Approx(20.0).margin(1e-9));
}

TEST_CASE("Heston pricer stays finite when the Feller condition is violated", "[heston][edge_case]") {
    // Real BTC vol surfaces routinely violate the Feller condition
    // (2*kappa*theta > xi^2, which would keep the variance process away
    // from zero); the COS method doesn't simulate the SDE path-by-path,
    // so it shouldn't care, but this is worth pinning down explicitly
    // since it's exactly the regime the live calibration will hit.
    deriv::HestonParams params{0.5, 0.3, 0.5, 2.5, -0.8};  // xi^2=6.25 >> 2*kappa*theta=0.3
    REQUIRE(2.0 * params.kappa * params.theta < params.xi * params.xi);

    deriv::MarketData market{100000.0, 0.05, 0.0, 0.0};
    deriv::EuropeanOption option{100000.0, 30.0 / 365.0, deriv::OptionType::Call};

    double price = deriv::heston_price(option, market, params);
    REQUIRE(std::isfinite(price));
    REQUIRE(price > 0.0);
}

TEST_CASE("Heston Monte Carlo (independent SDE simulation) agrees with the COS-method closed form",
          "[heston][monte_carlo][property]") {
    // heston.cpp's own header comment names two independent correctness
    // checks for the pricer: the deterministic-variance limit against
    // Black-Scholes (see the property test above), and an actual Monte
    // Carlo simulation of the Heston SDEs agreeing with the closed-form
    // COS price. This is that second check, finally implemented.
    //
    // It matters precisely because it's a DIFFERENT kind of check than
    // everything else in this file: every other Heston test exercises
    // heston_log_return_cf() and the COS summation -- if those share a
    // bug (e.g. a sign convention both heston_price() and a hypothetical
    // wrong derivation would agree on), no amount of testing heston.cpp
    // against itself would catch it. heston_monte_carlo_price() shares
    // NONE of that code: it simulates dS and dv directly via full-
    // truncation Euler (see heston_mc.hpp) and prices by discounted
    // average payoff. Agreement between the two is real evidence the CF
    // derivation matches the model it's supposed to be the CF of.
    //
    // Path/step counts (40,000 paths x 100 steps, antithetic) and the
    // fixed seed were chosen empirically: this combination keeps the
    // observed COS-vs-MC gap within about 1 standard error across many
    // seeds and scenarios (checked manually, not asserted here, since
    // that would just be re-deriving the CLT) while running in well
    // under a second per case. The margin below (5 standard errors) is
    // deliberately generous -- under 1e-6 false-failure probability for
    // a correctly-implemented pricer -- so this test fails only for a
    // real regression, not sampling noise.
    struct Scenario {
        double spot, strike, T;
        deriv::HestonParams heston;
        deriv::OptionType type;
        const char* label;
        int num_paths = 40000;
        int num_steps = 100;
    };

    std::vector<Scenario> scenarios = {
        {100.0, 100.0, 1.0, {0.04, 2.0, 0.04, 0.5, -0.5}, deriv::OptionType::Call, "ATM, 1yr, moderate vol-of-vol"},
        {100000.0, 120000.0, 90.0 / 365.0, {0.64, 1.5, 0.64, 0.9, -0.6}, deriv::OptionType::Call,
         "BTC-like, 90d, high vol-of-vol"},
        {100.0, 90.0, 2.0, {0.09, 3.0, 0.16, 0.3, 0.4}, deriv::OptionType::Put,
         "positive rho, v0 != theta, 2yr put"},
        {100.0, 100.0, 0.25, {0.25, 4.0, 0.25, 1.2, -0.7}, deriv::OptionType::Call,
         "short-dated, extreme vol-of-vol"},
        // The live demo's own deep-OTM test case (strike = 3x spot) with
        // its actual illustrative smile params (kappa=2, xi=0.5,
        // rho=-0.6): the ~$0.07 COS price here looked suspiciously far
        // below Black-Scholes' ~$1.65 at first glance, but this
        // independent SDE simulation confirms it's correct, not a bug --
        // negative spot-vol correlation genuinely thins the right tail
        // relative to constant-vol Black-Scholes, so a deep OTM call is
        // legitimately cheaper under Heston. See also the extreme-
        // moneyness case added to the deterministic-variance-limit test
        // above, which pins down a real (separate) precision bug this
        // same scenario exposed in double-precision put-call parity.
        //
        // This option finishes in-the-money on only a tiny fraction of
        // paths, so it needs far more than the default 40,000 paths to
        // get a nonzero (let alone informative) standard error -- 500,000
        // paths keeps this well under a second while reliably sampling
        // enough ITM paths to compare against the COS price.
        {100000.0, 300000.0, 90.0 / 365.0, {0.36, 2.0, 0.36, 0.5, -0.6}, deriv::OptionType::Call,
         "BTC deep OTM call, strike=3x spot, matches live demo smile params", 500000, 150},
    };

    for (const auto& s : scenarios) {
        INFO(s.label);

        deriv::MarketData market{s.spot, 0.05, 0.0, 0.0};
        deriv::EuropeanOption option{s.strike, s.T, s.type};

        double cos_price = deriv::heston_price(option, market, s.heston);
        deriv::MonteCarloResult mc =
            deriv::heston_monte_carlo_price(option, market, s.heston, s.num_paths, s.num_steps, true, 42);

        INFO("COS price: " << cos_price << ", MC price: " << mc.price << " +/- " << mc.standard_error);
        REQUIRE(mc.standard_error > 0.0);
        REQUIRE(cos_price == Catch::Approx(mc.price).margin(5.0 * mc.standard_error));
    }
}

// --- Heston calibration -----------------------------------------------

TEST_CASE("Heston calibration recovers a synthetic smile from a deliberately bad initial guess",
          "[heston][calibration]") {
    // Ground truth: generate a full synthetic options chain (12 strikes x
    // both call/put, one expiry) FROM a known Heston parameter set, solve
    // each contract's Black-Scholes-equivalent implied vol the same way
    // build_vol_surface_from_snapshot() would for real data, then check
    // that calibrate_heston() -- starting from parameters that don't
    // resemble the truth at all -- converges to a fit that reprices this
    // synthetic smile accurately.
    //
    // This deliberately does NOT assert the calibrated (v0, kappa, theta,
    // xi, rho) match the ground-truth values component-by-component.
    // That would be the wrong test: Heston's kappa (mean-reversion speed)
    // and theta (long-run variance) are not separately identifiable from
    // a single expiry's smile -- many (kappa, theta) pairs price that one
    // maturity almost identically, since what actually shows up in a
    // single-T price is closer to an integrated combination of the two,
    // not either one alone. Multiple *equally valid* parameter sets can
    // fit the same single-expiry smile; asserting one specific set is the
    // "right" answer would make this test fail on a correct calibrator.
    // What's actually well-identified from one expiry -- v0 (sets the
    // ATM level) and rho (sets the skew direction/steepness) -- IS
    // checked directly below, and the fit quality (which is the thing
    // that actually matters for repricing) is checked for everyone.
    deriv::HestonParams truth{0.36, 2.0, 0.30, 0.7, -0.55};
    double spot = 100000.0, r = 0.05, T = 90.0 / 365.0;

    std::vector<double> strikes = {70000, 80000, 85000, 90000, 95000, 100000,
                                    105000, 110000, 115000, 120000, 130000, 150000};
    std::vector<deriv::VolSurfacePoint> points;
    for (double K : strikes) {
        for (deriv::OptionType type : {deriv::OptionType::Call, deriv::OptionType::Put}) {
            deriv::EuropeanOption opt{K, T, type};
            deriv::MarketData market{spot, r, 0.0, 0.0};
            double price = deriv::heston_price(opt, market, truth);
            deriv::ImpliedVolResult iv = deriv::implied_volatility(opt, market, price);

            deriv::VolSurfacePoint p;
            p.strike = K;
            p.time_to_expiry = T;
            p.moneyness = K / spot;
            p.type = type;
            p.market_iv_reported = iv.vol;
            p.solved_iv = iv.vol;
            p.converged = iv.converged;
            p.underlying_price = spot;
            p.market_price_usd = price;
            p.risk_free_rate = r;
            points.push_back(p);
        }
    }

    // Deliberately far from the truth: half the ATM variance, half the
    // mean-reversion speed, less than half the vol-of-vol, weaker
    // correlation. If calibration only "works" when started near the
    // answer, that's not really calibration.
    deriv::HestonParams bad_guess{0.09, 1.0, 0.09, 0.3, -0.3};

    deriv::HestonCalibrationResult result = deriv::calibrate_heston(points, spot, r, bad_guess);

    INFO("calibrated v0=" << result.params.v0 << " kappa=" << result.params.kappa
                           << " theta=" << result.params.theta << " xi=" << result.params.xi
                           << " rho=" << result.params.rho);
    INFO("rmse_relative_price=" << result.rmse_relative_price << " rmse_iv_pp=" << result.rmse_iv_pp);

    REQUIRE(result.num_points_used == static_cast<int>(points.size()));
    REQUIRE(result.rmse_iv_pp < 0.05);           // sub-0.05-percentage-point fit to the smile
    REQUIRE(result.rmse_relative_price < 0.005);  // sub-0.5% price fit

    // v0 and rho ARE identifiable from a single expiry (see comment
    // above) and should come back close to the ground truth.
    REQUIRE(result.params.v0 == Catch::Approx(truth.v0).epsilon(0.05));
    REQUIRE(result.params.rho == Catch::Approx(truth.rho).epsilon(0.05));
}

TEST_CASE("Heston calibrates to the real Deribit BTC chain with a competitive fit",
          "[.manual][heston][calibration][real_data]") {
    // Calibrates to a single expiry slice of the real, committed Deribit
    // snapshot -- the nearest-to-90-day expiry with reasonable liquidity
    // -- restricted to a moneyness band of [0.7, 1.4]. That restriction
    // is a deliberate, documented choice, not a hidden filter to make the
    // numbers look better: very deep ITM/OTM contracts on the real chain
    // have tiny open interest and wide effective spreads (the mark_price
    // is thinly supported), and because the calibration objective is a
    // *relative* price error (see calibrate_heston's header comment),
    // those thin far-wing quotes get outsized weight in the objective
    // relative to how much genuine smile information they carry --
    // pulling the fit toward noise instead of signal. Real trading desks
    // apply the same kind of liquidity/moneyness filter before
    // calibrating for exactly this reason.
    auto all_points = deriv::build_vol_surface_from_snapshot("data/btc_chain_snapshot.json", 2026, 9, 17);
    REQUIRE(all_points.size() > 500);  // sanity check the snapshot loaded

    double target_T = 90.0 / 365.0;
    double chosen_key = -1.0;
    {
        std::vector<std::pair<double, int>> expiry_counts;
        for (const auto& p : all_points) {
            double key = std::round(p.time_to_expiry * 3650.0) / 3650.0;  // bucket to ~0.1 day
            bool found = false;
            for (auto& kv : expiry_counts) {
                if (std::abs(kv.first - key) < 1e-9) {
                    kv.second++;
                    found = true;
                    break;
                }
            }
            if (!found) expiry_counts.push_back({key, 1});
        }
        double best_dist = 1e18;
        for (const auto& kv : expiry_counts) {
            if (kv.second < 10) continue;  // need enough contracts to calibrate against
            double dist = std::abs(kv.first - target_T);
            if (dist < best_dist) {
                best_dist = dist;
                chosen_key = kv.first;
            }
        }
    }
    REQUIRE(chosen_key > 0.0);

    std::vector<deriv::VolSurfacePoint> slice;
    double spot = 0.0, r = 0.0;
    for (const auto& p : all_points) {
        double key = std::round(p.time_to_expiry * 3650.0) / 3650.0;
        if (std::abs(key - chosen_key) < 1e-9 && p.converged && p.moneyness > 0.7 && p.moneyness < 1.4) {
            slice.push_back(p);
            spot = p.underlying_price;
            r = p.risk_free_rate;
        }
    }
    REQUIRE(slice.size() > 20);

    double atm_iv = 0.6;
    double best_moneyness_dist = 1e18;
    for (const auto& p : slice) {
        double dist = std::abs(p.moneyness - 1.0);
        if (dist < best_moneyness_dist) {
            best_moneyness_dist = dist;
            atm_iv = p.solved_iv;
        }
    }

    deriv::HestonParams guess{atm_iv * atm_iv, 2.0, atm_iv * atm_iv, 0.7, -0.5};
    // regularization_weight=0.001 matches tools/calibrate_heston.cpp (the
    // program that actually feeds the live demo) -- see its comment, and
    // calibrate_heston's header comment, for why: kappa/theta aren't
    // separately identifiable from one expiry, and this keeps the fit
    // out of the implausible-but-technically-valid corner of parameter
    // space an unregularized fit was observed to find.
    deriv::HestonCalibrationResult result = deriv::calibrate_heston(slice, spot, r, guess, 0.001);

    WARN("Expiry: " << chosen_key * 365.0 << " days, " << slice.size() << " contracts in [0.7, 1.4] moneyness");
    WARN("Calibrated: v0=" << result.params.v0 << " kappa=" << result.params.kappa
                            << " theta=" << result.params.theta << " xi=" << result.params.xi
                            << " rho=" << result.params.rho);
    WARN("Fit quality: rmse_iv_pp=" << result.rmse_iv_pp
                                     << ", rmse_relative_price=" << result.rmse_relative_price);

    REQUIRE(result.converged);
    // A single flat vol (Phase 8's headline result) averages ~1.9
    // percentage points from Deribit's own reported IV across the whole
    // chain. Heston, with 5 free parameters fit to one expiry's actual
    // smile shape, should comfortably beat that on its own slice --
    // generously bounded at 1.0pp so this doesn't become flaky on a
    // future snapshot with a different real smile shape, while still
    // catching an actual calibration regression (a broken fit here lands
    // at several percentage points or fails to converge at all, not
    // hovers just above 1.0).
    REQUIRE(result.rmse_iv_pp < 1.0);
}

TEST_CASE("Multi-expiry Heston calibration resolves the single-expiry kappa/theta "
          "non-identifiability",
          "[heston][calibration][multi_expiry]") {
    // calibrate_heston() is generic over expiries -- see heston_calibration.hpp's
    // updated header comment -- but the single-expiry test above deliberately does
    // NOT check kappa/theta against ground truth, because one expiry's smile does
    // not pin them down uniquely. This test makes that identifiability gap, and
    // its resolution, directly visible: the SAME synthetic ground truth, the SAME
    // deliberately-bad initial guess, calibrated once against a single expiry and
    // once against three pooled expiries.
    //
    // v0 != theta is essential here: it creates a genuine variance term structure
    // (v(t) decays from v0 toward theta as t grows) for multiple expiries to carry
    // different information about. If v0 == theta, variance never mean-reverts and
    // kappa stays unidentifiable no matter how many expiries are used.
    deriv::HestonParams truth{0.16, 2.0, 0.09, 0.5, -0.6};
    double spot = 50000.0, r = 0.03;
    std::vector<double> strikes = {40000, 45000, 47500, 50000, 52500, 55000, 60000};
    std::vector<double> Ts = {30.0 / 365.0, 90.0 / 365.0, 270.0 / 365.0};

    auto gen_points = [&](double T) {
        std::vector<deriv::VolSurfacePoint> pts;
        for (double K : strikes) {
            for (deriv::OptionType type : {deriv::OptionType::Call, deriv::OptionType::Put}) {
                deriv::EuropeanOption opt{K, T, type};
                deriv::MarketData market{spot, r, 0.0, 0.0};
                double price = deriv::heston_price(opt, market, truth);
                deriv::ImpliedVolResult iv = deriv::implied_volatility(opt, market, price);

                deriv::VolSurfacePoint p;
                p.strike = K;
                p.time_to_expiry = T;
                p.moneyness = K / spot;
                p.type = type;
                p.market_iv_reported = iv.vol;
                p.solved_iv = iv.vol;
                p.converged = iv.converged;
                p.underlying_price = spot;
                p.market_price_usd = price;
                p.risk_free_rate = r;
                pts.push_back(p);
            }
        }
        return pts;
    };

    std::vector<deriv::VolSurfacePoint> single_expiry = gen_points(Ts[1]);  // 90d only
    std::vector<deriv::VolSurfacePoint> pooled;
    for (double T : Ts) {
        auto pts = gen_points(T);
        pooled.insert(pooled.end(), pts.begin(), pts.end());
    }

    // v0 correct, kappa >3x too fast, theta correct, xi and rho both off --
    // deliberately not "near the answer".
    deriv::HestonParams bad_guess{0.16, 6.0, 0.09, 0.3, -0.4};

    deriv::HestonCalibrationResult single = deriv::calibrate_heston(single_expiry, spot, r, bad_guess, 0.0);
    deriv::HestonCalibrationResult multi = deriv::calibrate_heston(pooled, spot, r, bad_guess, 0.0);

    INFO("single-expiry: v0=" << single.params.v0 << " kappa=" << single.params.kappa
                               << " theta=" << single.params.theta << " xi=" << single.params.xi
                               << " rho=" << single.params.rho << " rmse_rel_price=" << single.rmse_relative_price);
    INFO("multi-expiry:  v0=" << multi.params.v0 << " kappa=" << multi.params.kappa
                               << " theta=" << multi.params.theta << " xi=" << multi.params.xi
                               << " rho=" << multi.params.rho << " rmse_rel_price=" << multi.rmse_relative_price);

    // Both fits should reprice their own data almost exactly -- this is what
    // makes the single-expiry result a genuine identifiability failure rather
    // than "just a worse fit": it's an EXCELLENT fit with the wrong kappa/theta.
    REQUIRE(single.rmse_relative_price < 0.01);
    REQUIRE(multi.rmse_relative_price < 0.01);

    // Single-expiry: despite that excellent fit, kappa is NOT recovered --
    // the optimizer is free to run off along the unidentified ridge.
    double single_kappa_rel_err = std::abs(single.params.kappa - truth.kappa) / truth.kappa;
    REQUIRE(single_kappa_rel_err > 0.3);

    // Multi-expiry: the SAME bad initial guess, but pooling three maturities
    // pins kappa and theta down to close to the true values.
    REQUIRE(multi.params.kappa == Catch::Approx(truth.kappa).epsilon(0.10));
    REQUIRE(multi.params.theta == Catch::Approx(truth.theta).epsilon(0.10));
    REQUIRE(multi.params.v0 == Catch::Approx(truth.v0).epsilon(0.05));
    REQUIRE(multi.params.rho == Catch::Approx(truth.rho).epsilon(0.05));
}
TEST_CASE("Geometric Asian option: Monte Carlo path simulation matches the Kemna-Vorst closed form",
          "[exotic][asian][monte_carlo]") {
    // Arithmetic Asian options have no closed form, so this test validates
    // the path-simulation machinery itself (shared by the arithmetic
    // pricer) against a case that DOES have one: the geometric average.
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::AsianOption option{100.0, 1.0, deriv::OptionType::Call};
    int num_steps = 50;

    double closed_form = deriv::geometric_asian_price_closed_form(option, market, num_steps);
    deriv::AsianResult mc = deriv::geometric_asian_price_mc(option, market, 300000, num_steps, 42);

    INFO("closed-form: " << closed_form << "  MC: " << mc.price << " +- " << mc.standard_error);
    REQUIRE(mc.price == Catch::Approx(closed_form).margin(4.0 * mc.standard_error));
}

TEST_CASE("Arithmetic Asian call price is at least the geometric Asian call price (AM-GM)",
          "[exotic][asian][monte_carlo]") {
    // Pathwise, arithmetic average >= geometric average (AM-GM inequality),
    // and max(x-K,0) is nondecreasing in x, so this holds path-by-path and
    // therefore in expectation -- a model-independent inequality, not just
    // a numerical coincidence for this particular market.
    deriv::MarketData market{100.0, 0.05, 0.0, 0.3};
    deriv::AsianOption option{100.0, 1.0, deriv::OptionType::Call};
    int num_steps = 50;
    unsigned int seed = 123;

    deriv::AsianResult arithmetic = deriv::asian_option_price_mc(option, market, 200000, num_steps, seed);
    deriv::AsianResult geometric = deriv::geometric_asian_price_mc(option, market, 200000, num_steps, seed);

    INFO("arithmetic: " << arithmetic.price << "  geometric: " << geometric.price);
    REQUIRE(arithmetic.price > geometric.price);
}

TEST_CASE("Asian call price is below the vanilla European call price", "[exotic][asian][monte_carlo]") {
    // Averaging the underlying's path dampens its effective volatility
    // (the average of many draws has lower variance than a single terminal
    // draw), so an Asian call is systematically cheaper than the vanilla
    // European call with the same strike/expiry.
    deriv::MarketData market{100.0, 0.05, 0.0, 0.3};
    deriv::AsianOption asian{100.0, 1.0, deriv::OptionType::Call};
    deriv::EuropeanOption vanilla{100.0, 1.0, deriv::OptionType::Call};

    deriv::AsianResult asian_result = deriv::asian_option_price_mc(asian, market, 200000, 50, 7);
    double vanilla_price = deriv::black_scholes_price(vanilla, market);

    INFO("Asian: " << asian_result.price << "  vanilla: " << vanilla_price);
    REQUIRE(asian_result.price < vanilla_price);
}

TEST_CASE("Asian option with a single monitoring date matches vanilla European Monte Carlo",
          "[exotic][asian][monte_carlo]") {
    // With num_steps=1, the "average" is sampled at exactly one date (T),
    // so avg(S) == S_T and the Asian payoff collapses to the vanilla
    // European payoff exactly. Both pricers step the same GBM dynamics
    // with the same seed, so they should agree closely.
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::AsianOption asian{100.0, 1.0, deriv::OptionType::Call};
    deriv::EuropeanOption vanilla{100.0, 1.0, deriv::OptionType::Call};
    unsigned int seed = 99;

    deriv::AsianResult asian_result = deriv::asian_option_price_mc(asian, market, 200000, 1, seed);
    deriv::MonteCarloResult vanilla_result =
        deriv::monte_carlo_price(vanilla, market, 200000, /*antithetic=*/false, seed);

    INFO("Asian(1 step): " << asian_result.price << "  vanilla MC: " << vanilla_result.price);
    REQUIRE(asian_result.price ==
            Catch::Approx(vanilla_result.price).margin(4.0 * (asian_result.standard_error + vanilla_result.standard_error)));
}

TEST_CASE("Down-and-out barrier call: Monte Carlo matches the method-of-images closed form",
          "[exotic][barrier][monte_carlo]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::BarrierOption option{100.0, 80.0, 1.0, deriv::OptionType::Call, deriv::BarrierDirection::DownAndOut};

    double closed_form = deriv::down_and_out_call_price_closed_form(option, market);
    deriv::BarrierResult mc = deriv::barrier_option_price_mc(option, market, 300000, 100, 42);

    INFO("closed-form: " << closed_form << "  MC: " << mc.price << " +- " << mc.standard_error);
    REQUIRE(mc.price == Catch::Approx(closed_form).margin(4.0 * mc.standard_error));
}

TEST_CASE("In/out barrier parity: down-and-out + down-and-in == vanilla", "[exotic][barrier][monte_carlo]") {
    // Model-independent identity: a path either touches the barrier or it
    // doesn't, so a knock-out and its same-direction knock-in partition
    // every path exactly, and together they always reproduce the vanilla
    // payoff -- true regardless of the underlying's dynamics.
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption vanilla{100.0, 1.0, deriv::OptionType::Call};
    deriv::BarrierOption out_opt{100.0, 80.0, 1.0, deriv::OptionType::Call, deriv::BarrierDirection::DownAndOut};
    deriv::BarrierOption in_opt{100.0, 80.0, 1.0, deriv::OptionType::Call, deriv::BarrierDirection::DownAndIn};
    unsigned int seed = 55;

    deriv::BarrierResult out_result = deriv::barrier_option_price_mc(out_opt, market, 200000, 100, seed);
    deriv::BarrierResult in_result = deriv::barrier_option_price_mc(in_opt, market, 200000, 100, seed);
    double vanilla_price = deriv::black_scholes_price(vanilla, market);

    double combined_se = out_result.standard_error + in_result.standard_error;
    INFO("out: " << out_result.price << "  in: " << in_result.price
                  << "  sum: " << out_result.price + in_result.price << "  vanilla: " << vanilla_price);
    REQUIRE(out_result.price + in_result.price == Catch::Approx(vanilla_price).margin(4.0 * combined_se));
}

TEST_CASE("In/out barrier parity: up-and-out + up-and-in == vanilla", "[exotic][barrier][monte_carlo]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption vanilla{100.0, 1.0, deriv::OptionType::Call};
    deriv::BarrierOption out_opt{100.0, 120.0, 1.0, deriv::OptionType::Call, deriv::BarrierDirection::UpAndOut};
    deriv::BarrierOption in_opt{100.0, 120.0, 1.0, deriv::OptionType::Call, deriv::BarrierDirection::UpAndIn};
    unsigned int seed = 55;

    deriv::BarrierResult out_result = deriv::barrier_option_price_mc(out_opt, market, 200000, 100, seed);
    deriv::BarrierResult in_result = deriv::barrier_option_price_mc(in_opt, market, 200000, 100, seed);
    double vanilla_price = deriv::black_scholes_price(vanilla, market);

    double combined_se = out_result.standard_error + in_result.standard_error;
    INFO("out: " << out_result.price << "  in: " << in_result.price
                  << "  sum: " << out_result.price + in_result.price << "  vanilla: " << vanilla_price);
    REQUIRE(out_result.price + in_result.price == Catch::Approx(vanilla_price).margin(4.0 * combined_se));
}

// --- Bates model (Heston + Merton jumps) -----------------------------------
//
// bates.hpp's own header comment names three independent correctness
// checks, mirroring the pattern established for Heston: (1) jump_intensity
// = 0 must collapse EXACTLY to heston_price(), since the jump CF becomes
// exp(0) = 1 identically; (2) in the deterministic-variance limit (xi ->
// 0, v0 = theta), Bates must collapse to the closed-form Merton (1976)
// jump-diffusion price; (3) an independent Monte Carlo simulation of the
// actual Bates SDE + jump process must agree with the COS-method closed
// form. All three are implemented below.

namespace {

// Independent closed-form oracle for Merton (1976) jump-diffusion pricing,
// NOT calling into bates.cpp/heston.cpp at all: a Poisson-weighted sum of
// Black-Scholes prices, each evaluated at a jump-adjusted volatility and
// rate (the standard Merton series). n is capped at 50 terms, which for
// the lambda*T magnitudes used below makes the tail's contribution many
// orders of magnitude below double precision (Poisson tail past n=50 for
// lambda*T ~ O(1) is astronomically small).
//
// Subtlety this function specifically exists to get right (see
// docs/phase-log.md's Phase 14 entry for the full debugging story):
// black_scholes_price() called with rate=r_n uses r_n for BOTH generating
// the forward/terminal distribution AND discounting (e^{-r_n*T}). The
// correct Merton formula needs r_n only for the forward -- the true
// discount rate is always r. Naively calling black_scholes_price(rate=r_n)
// silently double-uses r_n and produces a price that is wrong by
// exp((r_n - r)*T), a small-looking (~0.5% in the case below) but real and
// systematic error that would make this "independent" oracle agree with
// nothing, correct or not. The fix: rescale each branch's raw
// black_scholes_price() output by exp((r_n - r)*T) before summing, which
// swaps the erroneous e^{-r_n*T} discount factor for the correct e^{-r*T}
// one while leaving the r_n-driven forward untouched.
double merton_jump_diffusion_price(double spot, double strike, double T, double r, double q, double vol,
                                    double lambda, double jump_mean, double jump_vol, deriv::OptionType type) {
    double k = std::exp(jump_mean + 0.5 * jump_vol * jump_vol) - 1.0;
    double sum = 0.0;

    for (int n = 0; n <= 50; ++n) {
        double sigma_n = std::sqrt(vol * vol + n * jump_vol * jump_vol / T);
        double r_n = r - lambda * k + n * std::log(1.0 + k) / T;

        deriv::MarketData market_n{spot, r_n, q, sigma_n};
        deriv::EuropeanOption option{strike, T, type};
        double bs_n = deriv::black_scholes_price(option, market_n);

        // Undo black_scholes_price()'s use of r_n as the discount rate;
        // see the function comment above.
        double bs_n_corrected = bs_n * std::exp((r_n - r) * T);

        double log_poisson_weight = -lambda * T + n * std::log(lambda * T) - std::lgamma(n + 1.0);
        sum += std::exp(log_poisson_weight) * bs_n_corrected;
    }

    return sum;
}

}  // namespace

TEST_CASE("Bates collapses exactly to Heston when jump intensity is zero", "[bates]") {
    // With jump_intensity = 0, jump_log_return_cf's compound-Poisson CF
    // is exp(0 * (phi - 1)) = exp(0) = 1 identically, and the compensator
    // k stops mattering (lambda*k = 0 regardless of jump_mean/jump_vol),
    // so bates_log_return_cf must reduce to heston_log_return_cf exactly
    // -- not approximately. This is the cheapest, most direct check that
    // the jump term was wired in correctly (multiplicatively, into the
    // CF) rather than, say, accidentally always contributing some residual
    // drift or variance.
    struct Scenario {
        double spot, strike, T;
        deriv::HestonParams heston;
        deriv::OptionType type;
        const char* label;
    };

    std::vector<Scenario> scenarios = {
        {100.0, 100.0, 1.0, {0.04, 2.0, 0.04, 0.5, -0.5}, deriv::OptionType::Call, "ATM call"},
        {100.0, 100.0, 1.0, {0.04, 2.0, 0.04, 0.5, -0.5}, deriv::OptionType::Put, "ATM put"},
        {100000.0, 120000.0, 90.0 / 365.0, {0.64, 1.5, 0.64, 0.9, -0.6}, deriv::OptionType::Call,
         "BTC-like OTM call"},
    };

    for (const auto& s : scenarios) {
        INFO(s.label);

        deriv::MarketData market{s.spot, 0.05, 0.0, 0.0};
        deriv::EuropeanOption option{s.strike, s.T, s.type};

        double heston = deriv::heston_price(option, market, s.heston);

        // Nonzero jump_mean/jump_vol deliberately, to confirm they're
        // inert when jump_intensity is zero -- not just that "all zeros"
        // happens to reduce correctly.
        deriv::BatesParams bates{s.heston, 0.0, -0.2, 0.4};
        double bates_p = deriv::bates_price(option, market, bates);

        REQUIRE(bates_p == Catch::Approx(heston).margin(1e-9));
    }
}

TEST_CASE("Bates collapses to the closed-form Merton jump-diffusion price in the "
          "deterministic-variance limit",
          "[bates][property]") {
    // Mirrors the Heston-vs-Black-Scholes deterministic-variance test
    // above: with xi -> 0 and v0 = theta, Heston's stochastic-vol
    // component degenerates to constant volatility, so Bates must
    // degenerate to plain Merton (1976) jump-diffusion -- priced here by
    // merton_jump_diffusion_price() above, which shares no code with
    // bates.cpp/heston.cpp.
    struct Scenario {
        double spot, strike, T, vol, lambda, jump_mean, jump_vol;
        deriv::OptionType type;
        const char* label;
    };

    std::vector<Scenario> scenarios = {
        {100.0, 100.0, 1.0, 0.20, 0.5, -0.1, 0.3, deriv::OptionType::Call, "ATM call, moderate jumps"},
        {100.0, 100.0, 1.0, 0.20, 0.5, -0.1, 0.3, deriv::OptionType::Put, "ATM put, moderate jumps"},
        {100.0, 120.0, 1.0, 0.20, 1.0, -0.15, 0.25, deriv::OptionType::Call, "OTM call, frequent jumps"},
        {100.0, 80.0, 0.5, 0.25, 0.3, 0.05, 0.2, deriv::OptionType::Put, "OTM put, positive mean jump, 6mo"},
        {100000.0, 110000.0, 90.0 / 365.0, 0.70, 2.0, -0.08, 0.35, deriv::OptionType::Call,
         "BTC-like: frequent jumps, 90d"},
    };

    // Same rationale as the Heston deterministic-variance test: xi small
    // but not literally zero, so the code path being tested still divides
    // by xi^2 (cancelled analytically) rather than skipping it.
    const double xi = 1e-6;
    const double r = 0.05, q = 0.0;

    for (const auto& s : scenarios) {
        INFO(s.label);

        deriv::HestonParams heston{s.vol * s.vol, 2.0, s.vol * s.vol, xi, -0.5};
        deriv::BatesParams bates{heston, s.lambda, s.jump_mean, s.jump_vol};

        deriv::MarketData market{s.spot, r, q, 0.0};
        deriv::EuropeanOption option{s.strike, s.T, s.type};

        double bates_p = deriv::bates_price(option, market, bates);
        double merton = merton_jump_diffusion_price(s.spot, s.strike, s.T, r, q, s.vol, s.lambda, s.jump_mean,
                                                      s.jump_vol, s.type);

        REQUIRE(bates_p == Catch::Approx(merton).epsilon(0.001));
    }
}

TEST_CASE("Bates Monte Carlo (independent SDE + jump simulation) agrees with the COS-method "
          "closed form",
          "[bates][monte_carlo][property]") {
    // The third leg of the validation triangle bates.hpp's header comment
    // names, mirroring heston_mc.hpp's role for Heston: bates_monte_carlo_price()
    // shares no code with bates_log_return_cf()/bates_put_price() -- it
    // simulates the Heston variance process (full-truncation Euler,
    // identical to heston_mc.cpp) with a compound Poisson jump added to
    // the log-price increment each step, and prices by discounted average
    // payoff. Agreement is real evidence the CF-based jump compensation
    // (r - lambda*k) and the compound-Poisson CF identity both match the
    // SDE they're supposed to describe.
    struct Scenario {
        double spot, strike, T;
        deriv::BatesParams bates;
        deriv::OptionType type;
        const char* label;
        int num_paths = 60000;
        int num_steps = 100;
    };

    std::vector<Scenario> scenarios = {
        {100.0, 100.0, 1.0, {{0.04, 2.0, 0.04, 0.5, -0.5}, 0.5, -0.1, 0.3}, deriv::OptionType::Call,
         "ATM, moderate vol-of-vol and jumps"},
        {100.0, 90.0, 2.0, {{0.09, 3.0, 0.16, 0.3, 0.4}, 0.3, 0.05, 0.2}, deriv::OptionType::Put,
         "positive rho, v0 != theta, mild positive-mean jumps, 2yr"},
        {100000.0, 120000.0, 90.0 / 365.0, {{0.64, 1.5, 0.64, 0.9, -0.6}, 1.0, -0.08, 0.35},
         deriv::OptionType::Call, "BTC-like: high vol-of-vol, frequent jumps, 90d"},
    };

    for (const auto& s : scenarios) {
        INFO(s.label);

        deriv::MarketData market{s.spot, 0.05, 0.0, 0.0};
        deriv::EuropeanOption option{s.strike, s.T, s.type};

        double cos_price = deriv::bates_price(option, market, s.bates);
        deriv::MonteCarloResult mc =
            deriv::bates_monte_carlo_price(option, market, s.bates, s.num_paths, s.num_steps, true, 42);

        INFO("COS price: " << cos_price << ", MC price: " << mc.price << " +/- " << mc.standard_error);
        REQUIRE(mc.standard_error > 0.0);
        REQUIRE(cos_price == Catch::Approx(mc.price).margin(5.0 * mc.standard_error));
    }
}

// --- Historical delta-hedging backtest (Phase 15) --------------------------

TEST_CASE("trailing_realized_vol matches a hand-computed value on a small series",
          "[backtest]") {
    // Three prices, two log-returns, checked against a value computed
    // independently by hand (not by re-running the same formula):
    // log(105/100) = 0.04879016417, log(103/105) = -0.01926909...
    // mean = 0.014760537, sample (n-1) variance over 2 points =
    // ((0.04879016417-0.014760537)^2 + (-0.01926909-0.014760537)^2) / 1
    // = (0.034029627)^2 + (-0.034029627)^2 = 2 * 0.034029627^2
    // daily_vol = sqrt(2 * 0.034029627^2) = 0.034029627 * sqrt(2)
    std::vector<deriv::PricePoint> history = {
        {"2020-01-01", 100.0}, {"2020-01-02", 105.0}, {"2020-01-03", 103.0}};

    double r0 = std::log(105.0 / 100.0);
    double r1 = std::log(103.0 / 105.0);
    double mean = (r0 + r1) / 2.0;
    double sq = (r0 - mean) * (r0 - mean) + (r1 - mean) * (r1 - mean);
    double sample_var = sq / 1.0;  // n-1 = 1
    double expected_daily_vol = std::sqrt(sample_var);
    double expected_annualized = expected_daily_vol * std::sqrt(365.0);

    double actual = deriv::trailing_realized_vol(history, 0, 2);
    REQUIRE(actual == Catch::Approx(expected_annualized).epsilon(1e-9));
}

TEST_CASE("trailing_realized_vol rejects an out-of-range window", "[backtest][edge_case]") {
    std::vector<deriv::PricePoint> history = {{"2020-01-01", 100.0}, {"2020-01-02", 101.0}};
    REQUIRE_THROWS_AS(deriv::trailing_realized_vol(history, 0, 5), std::runtime_error);
    REQUIRE_THROWS_AS(deriv::trailing_realized_vol(history, -1, 1), std::runtime_error);
}

namespace {

// Minimal Black-Scholes ModelFactory, mirroring tools/run_backtest.cpp's
// (not reused directly since that file's factories are local to its own
// anonymous namespace and this project's tools/ programs aren't part of
// the CMakeLists targets tests links against -- see tools/gen_stats.cpp
// and tools/calibrate_heston.cpp for the established precedent of
// standalone tools with their own main()). Kept intentionally tiny.
deriv::ModelFactory test_black_scholes_factory() {
    return [](double strike, double vol) -> deriv::PriceDeltaFn {
        return [strike, vol](double spot, double tau) -> std::pair<double, double> {
            deriv::EuropeanOption option{strike, tau, deriv::OptionType::Call};
            deriv::MarketData market{spot, 0.0, 0.0, vol};
            double price = deriv::black_scholes_price(option, market);
            deriv::Greeks g = deriv::black_scholes_greeks(option, market);
            return {price, g.delta};
        };
    };
}

deriv::ModelFactory test_heston_factory() {
    return [](double strike, double window_vol) -> deriv::PriceDeltaFn {
        double v0 = window_vol * window_vol;
        deriv::HestonParams params{v0, 2.0, v0, 0.5, -0.5};
        return [strike, params](double spot, double tau) -> std::pair<double, double> {
            deriv::EuropeanOption option{strike, tau, deriv::OptionType::Call};
            auto price_at = [&](double s) {
                deriv::MarketData m{s, 0.0, 0.0, 0.0};
                return deriv::heston_price(option, m, params);
            };
            double price = price_at(spot);
            double h = spot * 0.001;
            double delta = (price_at(spot + h) - price_at(spot - h)) / (2.0 * h);
            return {price, delta};
        };
    };
}

deriv::ModelFactory test_bates_factory(double jump_intensity) {
    return [jump_intensity](double strike, double window_vol) -> deriv::PriceDeltaFn {
        double v0 = window_vol * window_vol;
        deriv::BatesParams params{{v0, 2.0, v0, 0.5, -0.5}, jump_intensity, -0.1, 0.3};
        return [strike, params](double spot, double tau) -> std::pair<double, double> {
            deriv::EuropeanOption option{strike, tau, deriv::OptionType::Call};
            auto price_at = [&](double s) {
                deriv::MarketData m{s, 0.0, 0.0, 0.0};
                return deriv::bates_price(option, m, params);
            };
            double price = price_at(spot);
            double h = spot * 0.001;
            double delta = (price_at(spot + h) - price_at(spot - h)) / (2.0 * h);
            return {price, delta};
        };
    };
}

// Deterministic synthetic price series: a smooth-ish oscillation with a
// slight upward drift, entirely reproducible (no RNG), used to check the
// window-boundary bookkeeping (run_backtest's index arithmetic) rather
// than to make any claim about real market behavior -- the real-data
// test below is what does that.
std::vector<deriv::PricePoint> make_synthetic_series(int n, double s0 = 100.0) {
    std::vector<deriv::PricePoint> series;
    series.reserve(n);
    double S = s0;
    for (int i = 0; i < n; ++i) {
        char buf[16];
        std::snprintf(buf, sizeof(buf), "day-%04d", i);
        series.push_back(deriv::PricePoint{buf, S});
        S *= std::exp(0.0003 + 0.02 * std::sin(i * 0.3));
    }
    return series;
}

}  // namespace

TEST_CASE("run_backtest produces the expected number of non-overlapping windows",
          "[backtest]") {
    // 50 days of synthetic data, window_days=5, trailing_vol_days=5: the
    // first window starts at index 5 (once 5 days of trailing history
    // exist), and each subsequent window consumes 5 more days. The loop
    // condition (window_start_idx + window_days < n) with n=50 admits
    // window_start_idx = 5, 10, ..., 40 (44 + 5 = 49 < 50 is the last
    // one that fits) -- 8 windows.
    auto series = make_synthetic_series(50);
    auto summary = deriv::run_backtest(series, "test-BS", test_black_scholes_factory(), 5, 5);

    REQUIRE(summary.num_windows == 8);
    REQUIRE(summary.windows.size() == 8);

    // First window's bookkeeping: starts exactly at index 5 (trailing_vol_days),
    // strike is set at-the-money (== that day's spot).
    REQUIRE(summary.windows[0].date_start == series[5].date);
    REQUIRE(summary.windows[0].spot_start == Catch::Approx(series[5].price));
    REQUIRE(summary.windows[0].strike == Catch::Approx(series[5].price));

    // Second window starts exactly window_days later, not overlapping
    // the first.
    REQUIRE(summary.windows[1].date_start == series[10].date);
}

TEST_CASE("run_backtest throws when there isn't enough history for one window",
          "[backtest][edge_case]") {
    auto series = make_synthetic_series(10);
    REQUIRE_THROWS_AS(deriv::run_backtest(series, "test-BS", test_black_scholes_factory(), 30, 30),
                       std::runtime_error);
}

TEST_CASE("run_backtest stays finite and sane on a real market-scale synthetic series",
          "[backtest]") {
    auto series = make_synthetic_series(400, 30000.0);
    auto summary = deriv::run_backtest(series, "test-BS", test_black_scholes_factory(), 30, 30);

    REQUIRE(summary.num_windows > 0);
    REQUIRE(std::isfinite(summary.mean_pnl_bps));
    REQUIRE(std::isfinite(summary.stdev_pnl_bps));
    for (const auto& w : summary.windows) {
        REQUIRE(std::isfinite(w.hedge_pnl));
        REQUIRE(w.payoff >= 0.0);
        REQUIRE(w.theoretical_price >= 0.0);
    }
}

TEST_CASE("Backtest with Bates collapses to the Heston backtest when jump intensity is zero",
          "[backtest][bates]") {
    // The same exact-collapse property Bates itself has (see the
    // "[bates]"-tagged test above) should carry through the backtest
    // harness unchanged: with jump_intensity=0, test_bates_factory's
    // pricer is bit-for-bit the same computation as test_heston_factory's
    // (bates_price reduces to heston_price identically), so the two
    // backtests should agree on every window, not just on average.
    auto series = make_synthetic_series(200, 40000.0);

    auto heston_summary = deriv::run_backtest(series, "Heston", test_heston_factory(), 20, 20);
    auto bates_summary = deriv::run_backtest(series, "Bates (no jumps)", test_bates_factory(0.0), 20, 20);

    REQUIRE(heston_summary.num_windows == bates_summary.num_windows);
    for (size_t i = 0; i < heston_summary.windows.size(); ++i) {
        REQUIRE(bates_summary.windows[i].hedge_pnl == Catch::Approx(heston_summary.windows[i].hedge_pnl).margin(1e-6));
        REQUIRE(bates_summary.windows[i].theoretical_price ==
                Catch::Approx(heston_summary.windows[i].theoretical_price).margin(1e-9));
    }
}

TEST_CASE("Historical backtest on the real full BTC price history runs and produces finite, "
          "sane results for all three models",
          "[.manual][backtest][real_data]") {
    // Loads the real committed price history (data/btc_price_history.csv,
    // sourced from a public GitHub dataset -- see docs/numerics.md Phase
    // 15 for provenance and its honest limitations) and runs all three
    // models exactly as tools/run_backtest.cpp does, using the same real
    // calibrated Heston structural parameters
    // (tools/calibrate_heston.cpp's output against the committed Deribit
    // snapshot) and the same data-derived Bates jump parameters. Tagged
    // [.manual] like the other real-data tests in this file (Heston
    // calibration, implied vol) since it depends on a specific committed
    // data file's exact contents rather than being a self-contained
    // property check -- run explicitly with `./build/tests "[backtest]"`,
    // not part of the default suite.
    auto history = deriv::load_price_history("data/btc_price_history.csv");
    REQUIRE(history.size() > 1000);

    const double kappa = 1.829505, theta = 0.134505, xi = 1.280247, rho = -0.194560;
    const double jump_intensity = 3.4626461121463663;
    const double jump_mean = -0.000773771446447099;
    const double jump_vol = 0.27331845903521157;

    auto heston_factory = [=](double strike, double window_vol) -> deriv::PriceDeltaFn {
        double v0 = window_vol * window_vol;
        deriv::HestonParams params{v0, kappa, theta, xi, rho};
        return [strike, params](double spot, double tau) -> std::pair<double, double> {
            deriv::EuropeanOption option{strike, tau, deriv::OptionType::Call};
            auto price_at = [&](double s) {
                deriv::MarketData m{s, 0.0, 0.0, 0.0};
                return deriv::heston_price(option, m, params);
            };
            double price = price_at(spot);
            double h = spot * 0.001;
            double delta = (price_at(spot + h) - price_at(spot - h)) / (2.0 * h);
            return {price, delta};
        };
    };

    auto bates_factory = [=](double strike, double window_vol) -> deriv::PriceDeltaFn {
        double v0 = window_vol * window_vol;
        deriv::BatesParams params{{v0, kappa, theta, xi, rho}, jump_intensity, jump_mean, jump_vol};
        return [strike, params](double spot, double tau) -> std::pair<double, double> {
            deriv::EuropeanOption option{strike, tau, deriv::OptionType::Call};
            auto price_at = [&](double s) {
                deriv::MarketData m{s, 0.0, 0.0, 0.0};
                return deriv::bates_price(option, m, params);
            };
            double price = price_at(spot);
            double h = spot * 0.001;
            double delta = (price_at(spot + h) - price_at(spot - h)) / (2.0 * h);
            return {price, delta};
        };
    };

    auto bs_summary = deriv::run_backtest(history, "Black-Scholes", test_black_scholes_factory());
    auto heston_summary = deriv::run_backtest(history, "Heston", heston_factory);
    auto bates_summary = deriv::run_backtest(history, "Bates", bates_factory);

    INFO("Black-Scholes: mean=" << bs_summary.mean_pnl_bps << "bps stdev=" << bs_summary.stdev_pnl_bps
                                 << "bps over " << bs_summary.num_windows << " windows");
    INFO("Heston: mean=" << heston_summary.mean_pnl_bps << "bps stdev=" << heston_summary.stdev_pnl_bps
                          << "bps over " << heston_summary.num_windows << " windows");
    INFO("Bates: mean=" << bates_summary.mean_pnl_bps << "bps stdev=" << bates_summary.stdev_pnl_bps
                         << "bps over " << bates_summary.num_windows << " windows");

    REQUIRE(bs_summary.num_windows > 100);
    REQUIRE(std::isfinite(bs_summary.mean_pnl_bps));
    REQUIRE(std::isfinite(heston_summary.mean_pnl_bps));
    REQUIRE(std::isfinite(bates_summary.mean_pnl_bps));
}

// --- Portfolio, VaR, and stress testing (Phase 16) -------------------------

namespace {

deriv::PortfolioPricer test_bs_pricer() {
    return [](const deriv::EuropeanOption& opt, const deriv::MarketData& m) {
        return deriv::black_scholes_price(opt, m);
    };
}

}  // namespace

TEST_CASE("portfolio_value sums signed position values, options and underlying alike",
          "[portfolio]") {
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption call{100.0, 1.0, deriv::OptionType::Call};

    deriv::Portfolio book;
    book.positions.push_back(deriv::Position{deriv::Position::Kind::Option, call, 3.0, "long 3 calls"});
    book.positions.push_back(deriv::Position{deriv::Position::Kind::Option, call, -1.0, "short 1 call"});
    book.positions.push_back(
        deriv::Position{deriv::Position::Kind::Underlying, deriv::EuropeanOption{}, 2.0, "long 2 spot"});

    double call_price = deriv::black_scholes_price(call, market);
    double expected = (3.0 - 1.0) * call_price + 2.0 * market.spot;

    REQUIRE(deriv::portfolio_value(book, market, test_bs_pricer()) == Catch::Approx(expected).margin(1e-9));
}

TEST_CASE("portfolio_delta on a single-option portfolio matches Black-Scholes' own analytical delta",
          "[portfolio]") {
    // A one-position sanity check: portfolio_delta()'s finite-difference
    // approach should agree closely with the SAME option's own
    // closed-form delta -- confirming the finite-difference machinery
    // itself is correct before trusting it on Heston/Bates (which have
    // no closed-form delta to check against directly; see backtest.cpp's
    // identical choice for the same reason).
    deriv::MarketData market{100.0, 0.05, 0.0, 0.2};
    deriv::EuropeanOption call{110.0, 0.5, deriv::OptionType::Call};

    deriv::Portfolio book;
    book.positions.push_back(deriv::Position{deriv::Position::Kind::Option, call, 7.0, "7 calls"});

    deriv::Greeks g = deriv::black_scholes_greeks(call, market);
    double expected_delta = 7.0 * g.delta;

    REQUIRE(deriv::portfolio_delta(book, market, test_bs_pricer()) == Catch::Approx(expected_delta).epsilon(1e-4));
}

TEST_CASE("decay_time reduces every option's time-to-expiry and floors at zero, leaves underlying alone",
          "[portfolio][edge_case]") {
    deriv::Portfolio book;
    book.positions.push_back(deriv::Position{deriv::Position::Kind::Option,
                                               deriv::EuropeanOption{100.0, 0.1, deriv::OptionType::Call}, 1.0, "a"});
    book.positions.push_back(
        deriv::Position{deriv::Position::Kind::Underlying, deriv::EuropeanOption{}, 5.0, "spot"});

    // 0.1 years is ~36.5 days; decaying by 100 days should floor at 0,
    // not go negative.
    deriv::Portfolio decayed = deriv::decay_time(book, 100.0);

    REQUIRE(decayed.positions[0].option.time_to_expiry == Catch::Approx(0.0).margin(1e-12));
    REQUIRE(decayed.positions[1].quantity == Catch::Approx(5.0));  // underlying untouched
}

TEST_CASE("The example market-maker book is partially, not fully, delta-hedged",
          "[portfolio]") {
    // The whole point of this example (see portfolio.hpp's header
    // comment): a book with ~0 delta has ~0 first-order risk by
    // construction, which would make VaR/stress testing measure almost
    // nothing. Checked directly rather than assumed: the hedged book's
    // |delta| should be smaller than the unhedged book's (the spot
    // position is doing SOME hedging work) but still clearly nonzero
    // (it isn't doing ALL the hedging work).
    double spot = 50000.0;
    deriv::MarketData market{spot, 0.0, 0.0, 0.7};

    deriv::Portfolio book = deriv::make_example_market_maker_book(spot);
    deriv::Portfolio unhedged = book;
    unhedged.positions.pop_back();  // drop the spot hedge leg

    double delta_hedged = deriv::portfolio_delta(book, market, test_bs_pricer());
    double delta_unhedged = deriv::portfolio_delta(unhedged, market, test_bs_pricer());

    REQUIRE(std::abs(delta_unhedged) > 1e-6);            // the unhedged book has real delta risk
    REQUIRE(std::abs(delta_hedged) < std::abs(delta_unhedged));  // the hedge reduces it...
    REQUIRE(std::abs(delta_hedged) > 1e-6);               // ...without eliminating it
}

TEST_CASE("run_stress_test reproduces a real historical crash exactly, spot move and elapsed days alike",
          "[portfolio][real_data]") {
    auto history = deriv::load_price_history("data/btc_price_history.csv");
    double spot = 77000.0;
    deriv::MarketData market{spot, 0.0, 0.0, 0.7};
    deriv::Portfolio book = deriv::make_example_market_maker_book(spot);

    deriv::StressScenario covid{"2020 COVID crash", "2020-02-19", "2020-03-13"};
    deriv::StressResult r = deriv::run_stress_test(book, market, history, test_bs_pricer(), covid);

    // Hand-checked against the committed CSV directly (grep):
    // 2020-02-19,10189.9959829746 -> 2020-03-13,5800.2089048343
    REQUIRE(r.spot_start == Catch::Approx(10189.9959829746));
    REQUIRE(r.spot_end == Catch::Approx(5800.2089048343));
    REQUIRE(r.spot_pct_change == Catch::Approx(-0.4308).margin(0.001));
    // A real, large loss on a net-short-delta, short-gamma book facing a
    // ~43% crash -- sign and rough order of magnitude checked directly,
    // not just "finite".
    REQUIRE(r.pnl < 0.0);
}

TEST_CASE("run_stress_test throws on a date not present in the price history",
          "[portfolio][edge_case]") {
    auto history = deriv::load_price_history("data/btc_price_history.csv");
    double spot = 77000.0;
    deriv::MarketData market{spot, 0.0, 0.0, 0.7};
    deriv::Portfolio book = deriv::make_example_market_maker_book(spot);

    deriv::StressScenario bogus{"not a real date", "1999-01-01", "2020-03-13"};
    REQUIRE_THROWS_AS(deriv::run_stress_test(book, market, history, test_bs_pricer(), bogus), std::runtime_error);
}

TEST_CASE("Monte Carlo VaR is larger under Bates (with jumps) than under Heston (without) at 99% confidence",
          "[portfolio][var]") {
    // A jump component should make large moves MORE likely, not less --
    // so a jump-aware VaR should be at least as large as the equivalent
    // no-jump VaR on the same short-gamma book, especially at the 99%
    // tail where jumps matter most. This is the same "jumps add real
    // tail risk a diffusion-only model misses" story Phase 15's backtest
    // told, now checked directly on VaR instead of hedge P&L.
    double spot = 50000.0;
    deriv::MarketData market{spot, 0.0, 0.0, 0.7};
    deriv::Portfolio book = deriv::make_example_market_maker_book(spot);

    deriv::HestonParams heston{0.49, 1.8, 0.5, 1.0, -0.5};
    deriv::BatesParams bates_no_jumps{heston, 0.0, -0.1, 0.3};
    deriv::BatesParams bates_with_jumps{heston, 3.0, -0.1, 0.3};

    auto var_no_jumps = deriv::monte_carlo_var(book, market, bates_no_jumps, test_bs_pricer(), 20000, 1, 42);
    auto var_with_jumps = deriv::monte_carlo_var(book, market, bates_with_jumps, test_bs_pricer(), 20000, 1, 42);

    INFO("no jumps var_99=" << var_no_jumps.var_99 << " with jumps var_99=" << var_with_jumps.var_99);
    REQUIRE(var_with_jumps.var_99 >= var_no_jumps.var_99);
}

TEST_CASE("Historical VaR and Monte Carlo VaR broadly agree, and both clear the delta-normal "
          "baseline's known blind spot, on the real example book",
          "[.manual][portfolio][var][real_data]") {
    // The actual Phase 16 headline: real committed price history +
    // real calibrated Heston/Bates structural parameters (the same
    // ones Phase 14/15 used), on the shipped example book. Tagged
    // [.manual] like the other real-data-dependent tests in this file --
    // run explicitly with `./build/tests "[var]"`, not part of the
    // default suite, since it depends on a specific committed data
    // file's exact contents and takes noticeably longer than a unit
    // test (thousands of Monte Carlo paths x the full historical
    // series).
    auto history = deriv::load_price_history("data/btc_price_history.csv");
    double spot = 77000.0;
    deriv::MarketData market{spot, 0.0, 0.0, 0.75};
    deriv::Portfolio book = deriv::make_example_market_maker_book(spot);

    deriv::HestonParams heston{0.75 * 0.75, 1.829505, 0.134505, 1.280247, -0.194560};
    deriv::BatesParams bates{heston, 3.4626461121463663, -0.000773771446447099, 0.27331845903521157};

    auto hvar = deriv::historical_var(book, market, history, test_bs_pricer(), 1);
    auto mcvar = deriv::monte_carlo_var(book, market, bates, test_bs_pricer(), 20000, 1, 42);
    double dnvar_95 = deriv::delta_normal_var(book, market, test_bs_pricer(), 0.75, 1, 0.95);
    double dnvar_99 = deriv::delta_normal_var(book, market, test_bs_pricer(), 0.75, 1, 0.99);

    INFO("historical var95=" << hvar.var_95 << " var99=" << hvar.var_99);
    INFO("monte carlo var95=" << mcvar.var_95 << " var99=" << mcvar.var_99);
    INFO("delta-normal var95=" << dnvar_95 << " var99=" << dnvar_99);

    // Both repricing-based methods should be well within an order of
    // magnitude of each other -- a real, if loose, cross-validation
    // bound (not an exact-match assertion, since they're genuinely
    // different computations; see var.hpp).
    REQUIRE(hvar.var_99 / mcvar.var_99 < 10.0);
    REQUIRE(mcvar.var_99 / hvar.var_99 < 10.0);

    // The actual finding this phase exists to demonstrate: delta-normal
    // VaR, which only sees the book's linear (delta) risk, misses a
    // material fraction of the tail risk both repricing-based methods
    // (which see the book's real gamma) correctly price in.
    REQUIRE(dnvar_99 < hvar.var_99);
    REQUIRE(dnvar_99 < mcvar.var_99);
}
