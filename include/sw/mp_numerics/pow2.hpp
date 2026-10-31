#pragma once
// mp-numerics -- exact scaling by powers of two in any number type.

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>

namespace sw::mp_numerics {

/// x * 2^k in T, by factors of 2 small enough to be values of T. The partial
/// products move monotonically from x toward the result, so when both lie in
/// T's normal range every step is exact, even where 2^k itself, or 2^k beyond
/// double's exponent range, is not a value of T.
template <typename T>
T scale_pow2(T x, int k) {
    const int step = std::max(1, std::min(512, std::numeric_limits<T>::max_exponent - 2));
    while (k > step)  { x *= T(std::ldexp(1.0, step));  k -= step; }
    while (k < -step) { x *= T(std::ldexp(1.0, -step)); k += step; }
    return k == 0 ? x : T(x * T(std::ldexp(1.0, k)));
}

/// The binary exponent e of a positive finite x: 2^e <= x < 2^(e+1). For a
/// type other than the hardware floats it is read through double, so it may
/// be one too high where double(x) rounds up to a power of two, and it needs
/// x inside double's range; a wider type overloads it.
template <typename T>
int binary_exponent(const T& x) {
    if constexpr (std::is_floating_point_v<T>) return std::ilogb(x);
    else return std::ilogb(double(x));
}

/// 2^k in T, for k inside T's exponent range, even beyond double's.
template <typename T>
T pow2(int k) {
    return scale_pow2(T(1), k);
}

} // namespace sw::mp_numerics
