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
#include "deriv-engine/bates.hpp"
#include "deriv-engine/price_history.hpp"
#include "deriv-engine/backtest.hpp"
#include "deriv-engine/portfolio.hpp"
#include "deriv-engine/var.hpp"
#include "deriv-engine/pnl_attribution.hpp"
#include "deriv-engine/cva.hpp"

#include <emscripten/emscripten.h>
#include <vector>

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

// Real BTC price history, fetched and parsed once by the browser (see
// bridge_load_price_history below) and reused by every backtest/VaR/
// stress-test call after that -- a single global is fine here since the
// demo page only ever has one history loaded at a time, unlike the
// library code itself (backtest.hpp/var.hpp), which takes the history
// as a parameter and has no such assumption.
std::vector<PricePoint> g_history;

// The last backtest run's per-window results, kept around so the demo
// can draw a chart of hedge P&L across all 195 real windows without
// re-running the backtest or returning a huge array through the bridge
// in one call -- see bridge_backtest_window_count/pnl_bps below.
std::vector<BacktestWindowResult> g_last_backtest_windows;

// Every card in this demo prices the same way the Phase 16 test suite
// did (test_bs_pricer() in tests/test_main.cpp): Black-Scholes for
// portfolio_value/portfolio_delta, regardless of which model priced the
// option originally. See docs/numerics.md Phase 16 for why that's a
// reasonable choice for a portfolio-level risk number.
PortfolioPricer market_maker_pricer() {
    return [](const EuropeanOption& opt, const MarketData& m) { return black_scholes_price(opt, m); };
}

// Central finite-difference delta -- the same approach
// tools/run_backtest.cpp and backtest.cpp itself use for Heston/Bates,
// which have no closed-form Greeks in this engine (see backtest.hpp's
// header comment for why that's a deliberate scope decision).
double finite_diff_delta(double spot, const std::function<double(double)>& price_at_spot) {
    double h = spot * 0.001;
    return (price_at_spot(spot + h) - price_at_spot(spot - h)) / (2.0 * h);
}

ModelFactory make_bs_backtest_factory() {
    return [](double strike, double vol) -> PriceDeltaFn {
        return [strike, vol](double spot, double tau) -> std::pair<double, double> {
            EuropeanOption option{strike, tau, OptionType::Call};
            MarketData market{spot, 0.0, 0.0, vol};
            double price = black_scholes_price(option, market);
            Greeks g = black_scholes_greeks(option, market);
            return {price, g.delta};
        };
    };
}

ModelFactory make_heston_backtest_factory(double kappa, double theta, double xi, double rho) {
    return [=](double strike, double window_vol) -> PriceDeltaFn {
        double v0 = window_vol * window_vol;
        HestonParams params{v0, kappa, theta, xi, rho};
        return [strike, params](double spot, double tau) -> std::pair<double, double> {
            EuropeanOption option{strike, tau, OptionType::Call};
            auto price_at = [&](double s) {
                MarketData m{s, 0.0, 0.0, 0.0};
                return heston_price(option, m, params);
            };
            double price = price_at(spot);
            double delta = finite_diff_delta(spot, price_at);
            return {price, delta};
        };
    };
}

ModelFactory make_bates_backtest_factory(double kappa, double theta, double xi, double rho,
                                          double jump_intensity, double jump_mean, double jump_vol) {
    return [=](double strike, double window_vol) -> PriceDeltaFn {
        double v0 = window_vol * window_vol;
        BatesParams params{{v0, kappa, theta, xi, rho}, jump_intensity, jump_mean, jump_vol};
        return [strike, params](double spot, double tau) -> std::pair<double, double> {
            EuropeanOption option{strike, tau, OptionType::Call};
            auto price_at = [&](double s) {
                MarketData m{s, 0.0, 0.0, 0.0};
                return bates_price(option, m, params);
            };
            double price = price_at(spot);
            double delta = finite_diff_delta(spot, price_at);
            return {price, delta};
        };
    };
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

// ---------------------------------------------------------------------
// Bates (Phase 14): same COS method heston_price uses, extended with the
// compound Poisson jump term. v0/theta are variances, matching
// bridge_heston_price's convention above.
// ---------------------------------------------------------------------
EMSCRIPTEN_KEEPALIVE
double bridge_bates_price(double spot, double rate, double div,
                           double strike, double T, int type,
                           double v0, double kappa, double theta, double xi, double rho,
                           double jump_intensity, double jump_mean, double jump_vol) {
    BatesParams params{{v0, kappa, theta, xi, rho}, jump_intensity, jump_mean, jump_vol};
    return bates_price(make_euro(strike, T, type), make_market(spot, rate, div, 0.0), params);
}

// ---------------------------------------------------------------------
// Real BTC price history (Phase 15/16): the browser fetches
// docs/btc_price_history.csv as plain text and hands it here -- this
// project's own CSV parser (price_history.cpp) does the actual parsing,
// not a JS reimplementation. Returns the number of rows parsed, or -1 on
// a parse failure (malformed CSV) so the caller can show a real error
// instead of the module silently aborting.
// ---------------------------------------------------------------------
EMSCRIPTEN_KEEPALIVE
int bridge_load_price_history(const char* csv_text) {
    try {
        g_history = load_price_history_from_string(csv_text);
        g_last_backtest_windows.clear();
        return static_cast<int>(g_history.size());
    } catch (...) {
        g_history.clear();
        return -1;
    }
}

EMSCRIPTEN_KEEPALIVE
int bridge_history_size() {
    return static_cast<int>(g_history.size());
}

// ---------------------------------------------------------------------
// Historical delta-hedging backtest (Phase 15), run live against
// whatever g_history bridge_load_price_history() loaded. model: 0 =
// Black-Scholes, 1 = Heston, 2 = Bates (Heston/Bates params ignored for
// model 0). Writes [mean_pnl_bps, stdev_pnl_bps, num_windows] into out
// and returns 1 on success, 0 if the history isn't loaded yet or is too
// short for even one window -- the same "not enough history" case
// run_backtest() itself throws on, caught here rather than propagated.
// Also stashes the per-window results in g_last_backtest_windows for
// bridge_backtest_window_count/pnl_bps to draw a chart from.
// ---------------------------------------------------------------------
EMSCRIPTEN_KEEPALIVE
int bridge_run_backtest(int model, int window_days, int trailing_vol_days,
                         double kappa, double theta, double xi, double rho,
                         double jump_intensity, double jump_mean, double jump_vol,
                         double* out) {
    if (g_history.empty()) return 0;
    try {
        ModelFactory factory;
        std::string model_name;
        if (model == 0) {
            factory = make_bs_backtest_factory();
            model_name = "Black-Scholes";
        } else if (model == 1) {
            factory = make_heston_backtest_factory(kappa, theta, xi, rho);
            model_name = "Heston";
        } else {
            factory = make_bates_backtest_factory(kappa, theta, xi, rho, jump_intensity, jump_mean, jump_vol);
            model_name = "Bates";
        }

        BacktestSummary summary = run_backtest(g_history, model_name, factory, window_days, trailing_vol_days);
        out[0] = summary.mean_pnl_bps;
        out[1] = summary.stdev_pnl_bps;
        out[2] = static_cast<double>(summary.num_windows);
        g_last_backtest_windows = summary.windows;
        return 1;
    } catch (...) {
        return 0;
    }
}

EMSCRIPTEN_KEEPALIVE
int bridge_backtest_window_count() {
    return static_cast<int>(g_last_backtest_windows.size());
}

EMSCRIPTEN_KEEPALIVE
double bridge_backtest_window_pnl_bps(int index) {
    if (index < 0 || index >= static_cast<int>(g_last_backtest_windows.size())) return 0.0;
    return g_last_backtest_windows[index].hedge_pnl_bps;
}

// ---------------------------------------------------------------------
// Portfolio-level VaR (Phase 16) on the shipped example short-strangle
// book (make_example_market_maker_book), computed three independent
// ways: real historical simulation (needs g_history loaded), Monte
// Carlo under the Bates SDE, and the delta-normal baseline. Writes
// [hist_var95, hist_var99, hist_mean_pnl, mc_var95, mc_var99,
// mc_mean_pnl, delta_normal_var95, delta_normal_var99, unhedged_delta,
// hedged_delta] into out. Returns 1 on success, 0 if g_history isn't
// loaded yet.
// ---------------------------------------------------------------------
EMSCRIPTEN_KEEPALIVE
int bridge_run_var(double spot, double vol, int mc_paths, unsigned int seed,
                    double kappa, double theta, double xi, double rho,
                    double jump_intensity, double jump_mean, double jump_vol,
                    double* out) {
    if (g_history.empty()) return 0;
    try {
        MarketData market{spot, 0.0, 0.0, vol};
        Portfolio book = make_example_market_maker_book(spot);
        PortfolioPricer pricer = market_maker_pricer();

        VaRResult hvar = historical_var(book, market, g_history, pricer, 1);

        HestonParams heston{vol * vol, kappa, theta, xi, rho};
        BatesParams bates{heston, jump_intensity, jump_mean, jump_vol};
        VaRResult mcvar = monte_carlo_var(book, market, bates, pricer, mc_paths, 1, seed);

        double dn95 = delta_normal_var(book, market, pricer, vol, 1, 0.95);
        double dn99 = delta_normal_var(book, market, pricer, vol, 1, 0.99);

        Portfolio unhedged = book;
        unhedged.positions.pop_back();  // drop the spot hedge leg -- see portfolio.cpp

        out[0] = hvar.var_95;
        out[1] = hvar.var_99;
        out[2] = hvar.mean_pnl;
        out[3] = mcvar.var_95;
        out[4] = mcvar.var_99;
        out[5] = mcvar.mean_pnl;
        out[6] = dn95;
        out[7] = dn99;
        out[8] = portfolio_delta(unhedged, market, pricer);
        out[9] = portfolio_delta(book, market, pricer);
        return 1;
    } catch (...) {
        return 0;
    }
}

// ---------------------------------------------------------------------
// Real historical crash replay (Phase 16) on the example book. scenario:
// 0 = 2020 COVID crash, 1 = 2022 FTX collapse, 2 = the 2021-2022 bear
// market -- the same three, same order, as default_stress_scenarios()
// in var.cpp; the demo hardcodes their labels/dates in JS to match (see
// docs/index.html), the same way it already hardcodes barrier-direction
// labels to match instrument.hpp's BarrierDirection order. Writes
// [spot_start, spot_end, pct_change, pnl, value_start, value_end] into
// out. Returns 1 on success, 0 if g_history isn't loaded yet or the
// scenario's dates aren't in it.
// ---------------------------------------------------------------------
EMSCRIPTEN_KEEPALIVE
int bridge_run_stress_test(double spot, double vol, int scenario, double* out) {
    if (g_history.empty()) return 0;
    try {
        auto scenarios = default_stress_scenarios();
        if (scenario < 0 || scenario >= static_cast<int>(scenarios.size())) return 0;

        MarketData market{spot, 0.0, 0.0, vol};
        Portfolio book = make_example_market_maker_book(spot);
        StressResult r = run_stress_test(book, market, g_history, market_maker_pricer(), scenarios[scenario]);

        out[0] = r.spot_start;
        out[1] = r.spot_end;
        out[2] = r.spot_pct_change;
        out[3] = r.pnl;
        out[4] = r.portfolio_value_start;
        out[5] = r.portfolio_value_end;
        return 1;
    } catch (...) {
        return 0;
    }
}

// ---------------------------------------------------------------------
// P&L attribution (Phase 18), on the same short-gamma example book the
// VaR/stress-test card above uses (make_example_market_maker_book) --
// deliberately the SAME book, not a new one, so this card reads as "the
// desk explaining today's move on the book it's already looking at",
// not an unrelated toy. Doesn't need g_history loaded (attribute_pnl
// only needs two market snapshots, not a price series), so this card
// works independently of the backtest/VaR cards' history fetch.
// spot_move_pct is a fraction (0.08 = +8%), vol_change_pts is in vol
// POINTS (0.05 = +5 points, e.g. 75% -> 80%), matching attribute_pnl's
// own dVol convention. Writes [total_pnl, delta_pnl, gamma_pnl,
// vega_pnl, theta_pnl, unexplained_pnl, delta, gamma, vega, theta] into
// out. Always succeeds for finite inputs; try/catch kept for
// consistency with every other bridge function here.
// ---------------------------------------------------------------------
EMSCRIPTEN_KEEPALIVE
int bridge_run_pnl_attribution(double spot, double vol, double spot_move_pct, double vol_change_pts,
                                double days_elapsed, double* out) {
    try {
        MarketData before{spot, 0.0, 0.0, vol};
        MarketData after{spot * (1.0 + spot_move_pct), 0.0, 0.0, vol + vol_change_pts};
        Portfolio book = make_example_market_maker_book(spot);

        PnLAttributionResult r = attribute_pnl(book, before, after, market_maker_pricer(), days_elapsed);

        out[0] = r.total_pnl;
        out[1] = r.delta_pnl;
        out[2] = r.gamma_pnl;
        out[3] = r.vega_pnl;
        out[4] = r.theta_pnl;
        out[5] = r.unexplained_pnl;
        out[6] = r.delta;
        out[7] = r.gamma;
        out[8] = r.vega;
        out[9] = r.theta;
        return 1;
    } catch (...) {
        return 0;
    }
}

// ---------------------------------------------------------------------
// Basic CVA (Phase 18). Deliberately prices a single long 90-day
// at-the-money call "bought OTC from a risky counterparty" -- NOT the
// market-maker book above, which is mostly short options and would
// mostly demonstrate the trivial zero-exposure case (see cva.hpp's
// header comment: CVA only cares about POSITIVE exposure, what the
// counterparty owes YOU). hazard_rate/recovery_rate/horizon_years/
// num_time_steps/num_paths/seed and the Bates structural + jump
// parameters all forward straight into compute_cva(); num_time_steps is
// fixed at 12 (roughly monthly resolution) rather than exposed in the
// UI, to keep the card's inputs to what actually changes the answer in
// an interesting way. Writes [cva, expected_exposure_peak,
// expected_exposure_avg, position_value] into out. Returns 1 on
// success, 0 on invalid inputs (e.g. horizon_years <= 0).
// ---------------------------------------------------------------------
EMSCRIPTEN_KEEPALIVE
int bridge_run_cva(double spot, double vol, double hazard_rate, double recovery_rate, double horizon_years,
                    int num_paths, unsigned int seed, double kappa, double theta, double xi, double rho,
                    double jump_intensity, double jump_mean, double jump_vol, double* out) {
    try {
        MarketData market{spot, 0.0, 0.0, vol};

        EuropeanOption call{spot, 90.0 / 365.0, OptionType::Call};
        Portfolio book;
        book.positions.push_back(Position{Position::Kind::Option, call, 1.0, "long 1 90d ATM call (OTC)"});

        HestonParams heston{vol * vol, kappa, theta, xi, rho};
        BatesParams bates{heston, jump_intensity, jump_mean, jump_vol};

        const int num_time_steps = 12;
        CVAResult r = compute_cva(book, market, bates, market_maker_pricer(), hazard_rate, recovery_rate,
                                   horizon_years, num_time_steps, num_paths, seed);

        out[0] = r.cva;
        out[1] = r.expected_exposure_peak;
        out[2] = r.expected_exposure_avg;
        out[3] = portfolio_value(book, market, market_maker_pricer());
        return 1;
    } catch (...) {
        return 0;
    }
}

// Standard "credit triangle" approximation converting a CDS spread to a
// flat hazard rate -- see cva.hpp's header comment for what this is and
// isn't (a real desk bootstraps a full credit curve from CDS quotes
// across many tenors; this is the textbook single-point approximation).
// Exposed as its own bridge call, not computed in JS, so this
// conversion also comes out of the real compiled engine, not a
// JS-side reimplementation of the formula.
EMSCRIPTEN_KEEPALIVE
double bridge_hazard_rate_from_cds_spread(double cds_spread, double recovery_rate) {
    try {
        return hazard_rate_from_cds_spread(cds_spread, recovery_rate);
    } catch (...) {
        return -1.0;
    }
}

}  // extern "C"
