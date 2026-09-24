// Thin C-linkage bridge so Emscripten can export deriv-engine's real
// pricing functions to JavaScript with a stable, primitive-only ABI.
// This file adds NO new pricing logic — every call below forwards
// straight into the actual engine code in src/*.cpp.

#include "deriv-engine/market_data.hpp"
#include "deriv-engine/instrument.hpp"
#include "deriv-engine/black_scholes.hpp"
#include "deriv-engine/greeks.hpp"
#include "deriv-engine/tree_pricer.hpp"
#include "deriv-engine/monte_carlo.hpp"
#include "deriv-engine/implied_vol.hpp"
#include "deriv-engine/heston.hpp"
#include "deriv-engine/exotic_options.hpp"

#include <emscripten/emscripten.h>

using namespace deriv;

namespace {

MarketData make_market(double spot, double rate, double div, double vol) {
    MarketData m;
    m.spot = spot;
    m.risk_free_rate = rate;
    m.dividend_yield = div;
    m.volatility = vol;
    return m;
}

EuropeanOption make_euro(double strike, double T, int type) {
    EuropeanOption o;
    o.strike = strike;
    o.time_to_expiry = T;
    o.type = (type == 0) ? OptionType::Call : OptionType::Put;
    return o;
}

AmericanOption make_amer(double strike, double T, int type) {
    AmericanOption o;
    o.strike = strike;
    o.time_to_expiry = T;
    o.type = (type == 0) ? OptionType::Call : OptionType::Put;
    return o;
}

AsianOption make_asian(double strike, double T, int type) {
    AsianOption o;
    o.strike = strike;
    o.time_to_expiry = T;
    o.type = (type == 0) ? OptionType::Call : OptionType::Put;
    return o;
}

}  // namespace

extern "C" {

// Closed-form Black-Scholes price. type: 0 = call, 1 = put.
EMSCRIPTEN_KEEPALIVE
double bridge_bs_price(double spot, double rate, double div, double vol,
                        double strike, double T, int type) {
    return black_scholes_price(make_euro(strike, T, type),
                                make_market(spot, rate, div, vol));
}

// Analytical Greeks. Writes [delta, gamma, vega, theta, rho] into out.
EMSCRIPTEN_KEEPALIVE
void bridge_bs_greeks(double spot, double rate, double div, double vol,
                       double strike, double T, int type, double* out) {
    Greeks g = black_scholes_greeks(make_euro(strike, T, type),
                                     make_market(spot, rate, div, vol));
    out[0] = g.delta;
    out[1] = g.gamma;
    out[2] = g.vega;
    out[3] = g.theta;
    out[4] = g.rho;
}

// Binomial tree price, European exercise.
EMSCRIPTEN_KEEPALIVE
double bridge_tree_price_european(double spot, double rate, double div, double vol,
                                   double strike, double T, int type, int steps) {
    return binomial_tree_price(make_euro(strike, T, type),
                                make_market(spot, rate, div, vol), steps);
}

// Binomial tree price, American exercise (early-exercise aware).
EMSCRIPTEN_KEEPALIVE
double bridge_tree_price_american(double spot, double rate, double div, double vol,
                                   double strike, double T, int type, int steps) {
    return binomial_tree_price(make_amer(strike, T, type),
                                make_market(spot, rate, div, vol), steps);
}

// Trinomial tree price, European exercise.
EMSCRIPTEN_KEEPALIVE
double bridge_trinomial_price_european(double spot, double rate, double div, double vol,
                                        double strike, double T, int type, int steps) {
    return trinomial_tree_price(make_euro(strike, T, type),
                                 make_market(spot, rate, div, vol), steps);
}

// Trinomial tree price, American exercise (early-exercise aware).
EMSCRIPTEN_KEEPALIVE
double bridge_trinomial_price_american(double spot, double rate, double div, double vol,
                                        double strike, double T, int type, int steps) {
    return trinomial_tree_price(make_amer(strike, T, type),
                                 make_market(spot, rate, div, vol), steps);
}

// Monte Carlo price. Writes [price, standard_error] into out.
EMSCRIPTEN_KEEPALIVE
void bridge_mc_price(double spot, double rate, double div, double vol,
                      double strike, double T, int type, int num_paths,
                      int antithetic, unsigned int seed, double* out) {
    MonteCarloResult r = monte_carlo_price(make_euro(strike, T, type),
                                            make_market(spot, rate, div, vol),
                                            num_paths, antithetic != 0, seed);
    out[0] = r.price;
    out[1] = r.standard_error;
}

// Newton-Raphson implied vol solve. Writes [vol, converged(0/1), iterations] into out.
EMSCRIPTEN_KEEPALIVE
void bridge_implied_vol(double spot, double rate, double div,
                         double strike, double T, int type, double market_price,
                         double* out) {
    MarketData m = make_market(spot, rate, div, 0.0);  // vol unused by solver
    ImpliedVolResult r = implied_volatility(make_euro(strike, T, type), m, market_price);
    out[0] = r.vol;
    out[1] = r.converged ? 1.0 : 0.0;
    out[2] = static_cast<double>(r.iterations);
}

// Heston stochastic-vol price via the COS method. type: 0 = call, 1 = put.
// v0/theta are variances (vol^2), not vols -- callers (JS included) must
// square a vol before passing it in here, matching HestonParams' own
// convention (see heston.hpp).
EMSCRIPTEN_KEEPALIVE
double bridge_heston_price(double spot, double rate, double div,
                            double strike, double T, int type,
                            double v0, double kappa, double theta, double xi, double rho) {
    HestonParams params{v0, kappa, theta, xi, rho};
    // volatility is unused by heston_price -- Heston prices off params, not market.volatility.
    return heston_price(make_euro(strike, T, type),
                         make_market(spot, rate, div, 0.0), params);
}

// Asian option, Monte Carlo. geometric: 0 = arithmetic average (the real
// payoff), 1 = geometric average (validation-only -- see
// bridge_geometric_asian_closed_form below). Writes [price, standard_error].
EMSCRIPTEN_KEEPALIVE
void bridge_asian_price(double spot, double rate, double div, double vol,
                         double strike, double T, int type,
                         int num_paths, int num_steps, unsigned int seed,
                         int geometric, double* out) {
    AsianOption o = make_asian(strike, T, type);
    AsianResult r = (geometric != 0)
        ? geometric_asian_price_mc(o, make_market(spot, rate, div, vol), num_paths, num_steps, seed)
        : asian_option_price_mc(o, make_market(spot, rate, div, vol), num_paths, num_steps, seed);
    out[0] = r.price;
    out[1] = r.standard_error;
}

// Kemna-Vorst closed form for the geometric-average Asian option --
// exists to validate the Monte Carlo path simulation above, not as a
// price for the real (arithmetic-average) instrument.
EMSCRIPTEN_KEEPALIVE
double bridge_geometric_asian_closed_form(double spot, double rate, double div, double vol,
                                           double strike, double T, int type, int num_steps) {
    AsianOption o = make_asian(strike, T, type);
    return geometric_asian_price_closed_form(o, make_market(spot, rate, div, vol), num_steps);
}

// Barrier option, Monte Carlo. direction: 0=DownAndOut, 1=DownAndIn,
// 2=UpAndOut, 3=UpAndIn (matches BarrierDirection's declaration order in
// instrument.hpp). Writes [price, standard_error].
EMSCRIPTEN_KEEPALIVE
void bridge_barrier_price(double spot, double rate, double div, double vol,
                           double strike, double barrier, double T, int type,
                           int direction, int num_paths, int num_steps, unsigned int seed,
                           double* out) {
    BarrierOption o;
    o.strike = strike;
    o.barrier = barrier;
    o.time_to_expiry = T;
    o.type = (type == 0) ? OptionType::Call : OptionType::Put;
    o.direction = static_cast<BarrierDirection>(direction);
    BarrierResult r = barrier_option_price_mc(o, make_market(spot, rate, div, vol),
                                               num_paths, num_steps, seed);
    out[0] = r.price;
    out[1] = r.standard_error;
}

// Method-of-images closed form for a down-and-out CALL, valid only when
// barrier <= min(spot, strike) -- see docs/numerics.md Phase 13. Exists
// to validate bridge_barrier_price above, not as a general barrier formula.
EMSCRIPTEN_KEEPALIVE
double bridge_down_and_out_call_closed_form(double spot, double rate, double div, double vol,
                                             double strike, double barrier, double T) {
    BarrierOption o;
    o.strike = strike;
    o.barrier = barrier;
    o.time_to_expiry = T;
    o.type = OptionType::Call;
    o.direction = BarrierDirection::DownAndOut;
    return down_and_out_call_price_closed_form(o, make_market(spot, rate, div, vol));
}

}  // extern "C"
