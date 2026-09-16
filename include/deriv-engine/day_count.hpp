#pragma once

namespace deriv {

enum class DayCountConvention {
    Act365,
    Act360
};

double year_fraction(int days_between, DayCountConvention convention);

}  // namespace deriv