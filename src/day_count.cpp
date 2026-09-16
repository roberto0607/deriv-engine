#include "deriv-engine/day_count.hpp"

namespace deriv {

double year_fraction(int days_between, DayCountConvention convention) {
    const double denominator = (convention == DayCountConvention::Act365) ? 365.0 : 360.0;
    return static_cast<double>(days_between) / denominator;
}

}  // namespace deriv