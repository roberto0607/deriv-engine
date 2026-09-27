#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include "deriv-engine/monte_carlo.hpp"

using namespace deriv;

namespace {
// BTC-like inputs: spot, rate, dividend yield, vol
const MarketData mkt{65000.0, 0.04, 0.0, 0.55};
// At-the-money call, 3 months out
const EuropeanOption call{65000.0, 0.25, OptionType::Call};
}

TEST_CASE("Monte Carlo throughput", "[bench]") {
    BENCHMARK("MC price, 1M paths") {
        return monte_carlo_price(call, mkt, 1'000'000, false);
    };
}

TEST_CASE("AAD vs bump-and-reprice", "[bench]") {
    constexpr int N = 100'000;

    BENCHMARK("AAD: price + delta, vega, theta, rho") {
        return monte_carlo_greeks(call, mkt, N);
    };

    BENCHMARK("Bump: 1 base price + 4 bumped reprices") {
        const double h = 1e-4;
        MarketData s = mkt; s.spot *= (1 + h);
        MarketData v = mkt; v.volatility += h;
        MarketData r = mkt; r.risk_free_rate += h;
        EuropeanOption t = call; t.time_to_expiry -= h;

        double sum = monte_carlo_price(call, mkt, N, false, 42).price;
        sum += monte_carlo_price(call, s, N, false, 42).price;
        sum += monte_carlo_price(call, v, N, false, 42).price;
        sum += monte_carlo_price(call, r, N, false, 42).price;
        sum += monte_carlo_price(t, mkt, N, false, 42).price;
        return sum;
    };
}
