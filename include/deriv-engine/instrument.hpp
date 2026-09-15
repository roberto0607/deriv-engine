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

}  // namespace deriv