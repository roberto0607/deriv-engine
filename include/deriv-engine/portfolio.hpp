#pragma once

#include "market_data.hpp"
#include "instrument.hpp"
#include <functional>
#include <string>
#include <vector>

namespace deriv {

// A single position in a book: either an option (struck, dated, typed,
// same EuropeanOption every other pricer in this project already uses)
// or a position in the underlying itself (spot BTC). quantity is signed:
// positive is long, negative is short -- so "short 10 calls, struck at
// 1.1x spot" is quantity = -10 on an OptionType::Call position, and
// "hedged with 3 BTC of spot" is a separate Underlying position with
// quantity = 3.0. A single Portfolio is just a list of these; there's no
// separate "netting" step, since summing signed quantities in
// portfolio_value()/portfolio_delta() below does that implicitly.
struct Position {
    enum class Kind { Option, Underlying };

    Kind kind;
    EuropeanOption option;  // only meaningful when kind == Option
    double quantity;        // signed: + long, - short. For Underlying, units of spot.
    std::string label;      // for reporting only, e.g. "short 10x 30d 25-delta call"
};

struct Portfolio {
    std::vector<Position> positions;
};

// A model-agnostic pricer: given an option and the market it's priced
// under, returns its price. This is deliberately the SAME shape as
// backtest.hpp's PriceDeltaFn family (minus the delta, which
// portfolio_delta() below computes once for the whole portfolio via
// finite difference rather than needing each pricer to expose its own
// analytical Greeks) -- so a Black-Scholes, Heston, or Bates pricer can
// all be passed in unchanged, exactly as Phase 15's backtest did.
using PortfolioPricer = std::function<double(const EuropeanOption&, const MarketData&)>;

// Sums quantity * price(position) over every option position, plus
// quantity * market.spot over every underlying position. This is the
// portfolio's total mark-to-market value under `market` and `pricer` --
// the base case every VaR/stress-test computation in var.hpp compares a
// shocked scenario's value against.
double portfolio_value(const Portfolio& portfolio, const MarketData& market, const PortfolioPricer& pricer);

// Portfolio-level delta via a single central finite difference on the
// WHOLE portfolio's value (bump market.spot by 0.1%, reprice everything,
// divide) -- not a sum of each position's own analytical delta. This is
// deliberately the same choice backtest.cpp made for Heston/Bates (see
// its header comment): it works identically regardless of which pricer
// is plugged in, so a mixed portfolio priced partly under Black-Scholes
// and partly under Heston (via a PortfolioPricer that dispatches per
// position -- see the example book in portfolio.cpp) still gets one
// consistent delta number, without needing per-model analytical Greeks
// wired in everywhere.
double portfolio_delta(const Portfolio& portfolio, const MarketData& market, const PortfolioPricer& pricer);

// A realistic market-maker book, shipped as a worked example alongside
// the generic API above (the same pattern exotic_options.hpp used --
// worked examples in the live demo, not just an API with no instance of
// its own use). Short strangles/straddles across a few strikes and two
// expiries (collecting premium, the position a market maker who'd sold
// options to clients would actually be in), partially offset with a
// long spot position -- NOT fully delta-hedged, so the book has real,
// nonzero gamma risk left over, which is exactly the risk VaR and stress
// testing exist to measure. A fully delta-hedged book has ~0 risk by
// construction and wouldn't make an interesting VaR example.
Portfolio make_example_market_maker_book(double spot);

// Returns a copy of `portfolio` with every option position's
// time-to-expiry reduced by `calendar_days` days (floored at 0 --
// options don't go negative-dated, they just expire worthless/intrinsic
// once repriced at tau=0, which every pricer in this project already
// handles). Underlying positions are unaffected. Used by var.hpp's
// horizon-based P&L computations (historical VaR, Monte Carlo VaR, and
// the stress tests) to reflect that time actually passes over the
// horizon being measured, not just spot.
Portfolio decay_time(const Portfolio& portfolio, double calendar_days);

}  // namespace deriv
