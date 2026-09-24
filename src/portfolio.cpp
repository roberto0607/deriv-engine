#include "deriv-engine/portfolio.hpp"
#include <cmath>

namespace deriv {

double portfolio_value(const Portfolio& portfolio, const MarketData& market, const PortfolioPricer& pricer) {
    double total = 0.0;
    for (const auto& p : portfolio.positions) {
        if (p.kind == Position::Kind::Option) {
            total += p.quantity * pricer(p.option, market);
        } else {
            total += p.quantity * market.spot;
        }
    }
    return total;
}

double portfolio_delta(const Portfolio& portfolio, const MarketData& market, const PortfolioPricer& pricer) {
    double h = market.spot * 0.001;

    MarketData up = market;
    up.spot = market.spot + h;
    MarketData down = market;
    down.spot = market.spot - h;

    double v_up = portfolio_value(portfolio, up, pricer);
    double v_down = portfolio_value(portfolio, down, pricer);

    return (v_up - v_down) / (2.0 * h);
}

Portfolio make_example_market_maker_book(double spot) {
    Portfolio book;

    const double T1 = 30.0 / 365.0;  // near-dated leg
    const double T2 = 60.0 / 365.0;  // far-dated leg

    // Short 30-day at-the-money straddle x5: the classic "sold clients a
    // bunch of ATM vol" position -- maximum short gamma exactly where
    // spot is right now, which is precisely why it's the riskiest part
    // of this book.
    book.positions.push_back(
        Position{Position::Kind::Option, EuropeanOption{spot, T1, OptionType::Call}, -5.0, "short 5x 30d ATM call"});
    book.positions.push_back(
        Position{Position::Kind::Option, EuropeanOption{spot, T1, OptionType::Put}, -5.0, "short 5x 30d ATM put"});

    // Short 60-day 10%-OTM strangle x5: collects less premium per
    // contract than the straddle (further from the money) but adds
    // short gamma/vega further out the curve and away from spot's
    // current level, the way a real book accumulates risk from more
    // than one client trade.
    book.positions.push_back(Position{Position::Kind::Option, EuropeanOption{spot * 1.10, T2, OptionType::Call},
                                       -5.0, "short 5x 60d 110% strike call"});
    book.positions.push_back(Position{Position::Kind::Option, EuropeanOption{spot * 0.90, T2, OptionType::Put},
                                       -5.0, "short 5x 60d 90% strike put"});

    // Partial delta hedge: long spot, sized to roughly halve (not
    // eliminate) the book's short-option delta exposure. The unhedged
    // straddle+strangle combination above nets to about -1.06 BTC of
    // delta (computed via portfolio_delta(), not hand-derived -- see
    // docs/numerics.md Phase 16), so +0.5 BTC brings the book to roughly
    // -0.56 BTC of delta: partially, not fully, hedged. A fully-hedged
    // book has ~0 delta risk by construction and would make VaR/stress
    // testing measure almost nothing.
    book.positions.push_back(Position{Position::Kind::Underlying, EuropeanOption{}, 0.5, "partial spot hedge"});

    return book;
}

Portfolio decay_time(const Portfolio& portfolio, double calendar_days) {
    Portfolio decayed = portfolio;
    double dt = calendar_days / 365.0;
    for (auto& p : decayed.positions) {
        if (p.kind == Position::Kind::Option) {
            p.option.time_to_expiry = std::max(p.option.time_to_expiry - dt, 0.0);
        }
    }
    return decayed;
}

}  // namespace deriv
