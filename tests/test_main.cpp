#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "deriv-engine/black_scholes.hpp"

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