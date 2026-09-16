#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "deriv-engine/black_scholes.hpp"
#include "deriv-engine/greeks.hpp"

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