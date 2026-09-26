#include "deriv-engine/pnl_attribution.hpp"
#include <cmath>

namespace deriv {

double portfolio_gamma(const Portfolio& portfolio, const MarketData& market, const PortfolioPricer& pricer) {
    double h = market.spot * 0.001;

    MarketData up = market;
    up.spot = market.spot + h;
    MarketData down = market;
    down.spot = market.spot - h;

    double v_up = portfolio_value(portfolio, up, pricer);
    double v_mid = portfolio_value(portfolio, market, pricer);
    double v_down = portfolio_value(portfolio, down, pricer);

    return (v_up - 2.0 * v_mid + v_down) / (h * h);
}

double portfolio_vega(const Portfolio& portfolio, const MarketData& market, const PortfolioPricer& pricer) {
    double h = 0.001;  // 0.1 vol point, e.g. 60.0% -> 60.1%

    MarketData up = market;
    up.volatility = market.volatility + h;
    MarketData down = market;
    down.volatility = market.volatility - h;

    double v_up = portfolio_value(portfolio, up, pricer);
    double v_down = portfolio_value(portfolio, down, pricer);

    return (v_up - v_down) / (2.0 * h);
}

double portfolio_theta(const Portfolio& portfolio, const MarketData& market, const PortfolioPricer& pricer) {
    const double bump_days = 1.0;

    double v_now = portfolio_value(portfolio, market, pricer);
    Portfolio decayed = decay_time(portfolio, bump_days);
    double v_later = portfolio_value(decayed, market, pricer);

    double dt_years = bump_days / 365.0;
    // decay_time() moves calendar time FORWARD (time_to_expiry down), so
    // (v_later - v_now)/dt is directly d(price)/d(calendar time) = theta
    // under the sign convention documented in pnl_attribution.hpp --
    // no extra negation needed here, unlike a naive d(price)/d(T) read.
    return (v_later - v_now) / dt_years;
}

PnLAttributionResult attribute_pnl(const Portfolio& portfolio, const MarketData& market_before,
                                    const MarketData& market_after, const PortfolioPricer& pricer,
                                    double time_elapsed_days) {
    double value_before = portfolio_value(portfolio, market_before, pricer);

    Portfolio decayed = decay_time(portfolio, time_elapsed_days);
    double value_after = portfolio_value(decayed, market_after, pricer);

    double total_pnl = value_after - value_before;

    double delta = portfolio_delta(portfolio, market_before, pricer);
    double gamma = portfolio_gamma(portfolio, market_before, pricer);
    double vega = portfolio_vega(portfolio, market_before, pricer);
    double theta = portfolio_theta(portfolio, market_before, pricer);

    double dS = market_after.spot - market_before.spot;
    double dVol = market_after.volatility - market_before.volatility;
    double dt_years = time_elapsed_days / 365.0;

    double delta_pnl = delta * dS;
    double gamma_pnl = 0.5 * gamma * dS * dS;
    double vega_pnl = vega * dVol;
    double theta_pnl = theta * dt_years;

    double explained = delta_pnl + gamma_pnl + vega_pnl + theta_pnl;
    double unexplained_pnl = total_pnl - explained;

    return PnLAttributionResult{total_pnl, delta_pnl, gamma_pnl, vega_pnl, theta_pnl, unexplained_pnl,
                                 delta,     gamma,     vega,      theta};
}

}  // namespace deriv
