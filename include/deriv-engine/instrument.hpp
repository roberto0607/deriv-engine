#pragma once

namespace deriv {

enum class OptionType {
    Call,
    Put
};

struct EuropeanOption {
    double strike;
    double time_to_expiry;
    OptionType type;
};

struct AmericanOption {
    double strike;
    double time_to_expiry;
    OptionType type;
};

// Average-price (arithmetic or geometric, selected by which pricing
// function is called -- see exotic_options.hpp) Asian option. Payoff is
// max(avg(S) - K, 0) for a call, max(K - avg(S), 0) for a put, where
// avg(S) is the average of the underlying's price sampled at
// num_steps equally-spaced monitoring dates from just after today
// through expiry (the monitoring frequency is a pricing-function
// parameter, not part of the instrument itself, since it's a
// discretization choice, not a contractual term).
struct AsianOption {
    double strike;
    double time_to_expiry;
    OptionType type;
};

// Which side of the barrier knocks the option in or out, and which
// direction the barrier sits relative to spot. "Out" options start
// alive and are extinguished (permanently, to a worthless payoff) the
// first time the path touches the barrier; "In" options start dead
// and only become a live vanilla option if the path touches the
// barrier at some point before expiry. See exotic_options.hpp for the
// in/out parity identity this taxonomy makes possible to test:
// (same-direction out) + (same-direction in) == vanilla, always,
// model-independently.
enum class BarrierDirection {
    DownAndOut,
    DownAndIn,
    UpAndOut,
    UpAndIn
};

struct BarrierOption {
    double strike;
    double barrier;
    double time_to_expiry;
    OptionType type;
    BarrierDirection direction;
};

}  // namespace deriv