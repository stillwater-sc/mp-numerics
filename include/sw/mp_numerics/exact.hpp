#pragma once
// mp-numerics -- exact values of binary floating-point numbers.
//
// Accuracy studies need the error of a computed result measured without
// rounding the measurement: c^2 + s^2 - 1 for a posit<64,2> rotation is about
// 1e-18, below double's resolution of 1 + 1e-18. This header turns a value of
// any binary number type into a Shewchuk expansion (Universal's ereal), in which
// sums and products are exact, so such residuals can be formed exactly and
// rounded once, to double, for reporting.
//
// A value is read out by peeling: h = double(x) is x rounded to 53 bits, and
// x - T(h) is exact because both lie in the same binade and the difference
// needs fewer bits than x has. Repeating until nothing is left gives a
// non-overlapping sum of doubles equal to x. That needs T to be a binary format
// (cfloat, posit, linear takum, fixpnt, dd, qd, the hardware floats); a
// logarithmic type such as lns is not, and is converted through double instead.
//
// A value outside double's exponent range is first scaled by a power of two
// inside T, which is exact while the result stays in T's normal range. Callers
// pass the binary exponent they scaled their data by, and every quantity
// derived from the same data is scaled alike.

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>

#include <universal/number/ereal/ereal.hpp>
#include <universal/number/lns/lns_fwd.hpp>

#include <sw/mp_numerics/pow2.hpp>

namespace sw::mp_numerics {

/// An exact real: a Shewchuk expansion of up to 19 doubles. 19 is the most
/// limbs whose last one stays a normal double, about 1000 bits of span.
using exact = sw::universal::ereal<19>;

/// Whether T's values are dyadic rationals that peeling reads exactly. A
/// logarithmic type specializes this to false.
template <typename T>
struct is_binary_format : std::true_type {};

/// lns values are 2^(m / 2^rbits): read through double, which rounds them
/// far below the type's own precision.
template <unsigned nbits, unsigned rbits, typename bt, auto... xtra>
struct is_binary_format<sw::universal::lns<nbits, rbits, bt, xtra...>> : std::false_type {};

/// Thrown by to_exact for a value it cannot read: one that, scaled by 2^-k,
/// lies outside double's range, or (were T not binary) has no finite peeling.
struct unreadable_value : std::runtime_error {
    using std::runtime_error::runtime_error;
};

/// x * 2^-k formed in T, then converted to double. A NaN or an infinity in x
/// stays one, so this doubles as a finiteness test for scaled data.
template <typename T>
double scaled_double(const T& x, int k = 0) {
    return double(scale_pow2(x, -k));
}

/// The exact value of x * 2^-k. x must be finite, and x * 2^-k must land in
/// T's normal range for the scaling to be exact; both are the caller's to
/// ensure. Throws unreadable_value if peeling does not terminate: x * 2^-k is
/// beyond double's range, or below it with bits that no double holds.
template <typename T>
exact to_exact(const T& x, int k = 0) {
    T y = scale_pow2(x, -k);
    if constexpr (!is_binary_format<T>::value) {
        return exact(double(y));
    } else {
        exact acc(0.0);
        for (int i = 0; i < 64 && y != T(0); ++i) {
            const double h = double(y);
            acc += exact(h);
            y = y - T(h);
        }
        if (y != T(0))
            throw unreadable_value("to_exact: peeling did not terminate for x * 2^-k with double(x * 2^-k) = " +
                                     std::to_string(double(scale_pow2(x, -k))) + ", k = " + std::to_string(k));
        return acc;
    }
}

/// |e| rounded to double.
inline double abs_double(const exact& e) { return std::abs(double(e)); }

} // namespace sw::mp_numerics
