#pragma once
// mp-numerics -- exact scaling by powers of two in any number type.

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

#include <universal/number/dd/dd_fwd.hpp>
#include <universal/number/qd/qd_fwd.hpp>

namespace sw::mp_numerics {

/// x * 2^k in T, by factors of 2 small enough to be values of T. The partial
/// products move monotonically from x toward the result, so when both lie in
/// T's normal range every step is exact, even where 2^k itself, or 2^k beyond
/// double's exponent range, is not a value of T. A type whose range passes
/// double's first takes steps of 2^(512 * 2^i), formed by squaring 2^512 and
/// so exact, to keep the cost logarithmic in |k|.
///
/// dd and qd scale limb by limb with their own ldexp instead. Their product
/// is not safe near double's overflow threshold: Universal's Dekker split
/// rounds the high half of a value near DBL_MAX up to 2^1024, so
/// dd(DBL_MAX) * 2^-512 is a NaN.
template <typename T>
T scale_pow2(T x, int k) {
    if constexpr (std::is_same_v<T, sw::universal::dd> || std::is_same_v<T, sw::universal::qd>)
        return ldexp(x, k);
    constexpr int cap = std::numeric_limits<T>::max_exponent - 2;
    if constexpr (cap > 1024) {
        struct ladder { std::vector<T> up, down; std::vector<int> size; };
        static const ladder L = [] {
            ladder l;
            T u = T(0x1.0p512), d = T(0x1.0p-512);
            for (int n = 512; n <= cap; n *= 2) {
                l.up.push_back(u);
                l.down.push_back(d);
                l.size.push_back(n);
                if (n > cap / 2) break;
                u *= u;
                d *= d;
            }
            return l;
        }();
        for (std::size_t i = L.size.size(); i-- > 0;) {
            while (k >= L.size[i])  { x *= L.up[i];   k -= L.size[i]; }
            while (k <= -L.size[i]) { x *= L.down[i]; k += L.size[i]; }
        }
    }
    const int step = std::max(1, std::min(512, cap));
    while (k > step)  { x *= T(std::ldexp(1.0, step));  k -= step; }
    while (k < -step) { x *= T(std::ldexp(1.0, -step)); k += step; }
    return k == 0 ? x : T(x * T(std::ldexp(1.0, k)));
}

/// The binary exponent e of a positive finite x: 2^e <= x < 2^(e+1). For a
/// type other than the hardware floats it is read through double: x is first
/// brought inside double's normal range by exact power-of-two steps, and e is
/// lowered by one where double(x) rounded up to a power of two. So it holds
/// for a type whose range is wider than double's (cfloat<128>, cfloat<256>),
/// and for one with more precision (posit<64,2>, dd, qd).
template <typename T>
int binary_exponent(const T& x) {
    if constexpr (std::is_floating_point_v<T>) {
        return std::ilogb(x);
    } else {
        // The steps double up to 2^(2^20), so 64 iterations reach past any
        // type's range; more means x was not finite.
        if (!(x > T(0))) throw std::domain_error("binary_exponent: x must be positive");
        T y = x;
        int shift = 0, iterations = 0;
        auto guard = [&] {
            if (++iterations > 64) throw std::domain_error("binary_exponent: x must be finite");
        };
        // Steps that double while x stays out of range, and halve once a
        // step would carry it past double's range on the other side.
        for (int step = 512; !(double(y) < 0x1.0p1000); guard()) {
            const T z = scale_pow2(y, -step);
            if (step > 512 && double(z) < 0x1.0p-1000) { step /= 2; continue; }
            y = z; shift += step; step = std::min(2 * step, 1 << 20);
        }
        for (int step = 512; double(y) < 0x1.0p-1000; guard()) {
            const T z = scale_pow2(y, step);
            if (step > 512 && !(double(z) < 0x1.0p1000)) { step /= 2; continue; }
            y = z; shift -= step; step = std::min(2 * step, 1 << 20);
        }
        int e = std::ilogb(double(y));
        if (y < T(std::ldexp(1.0, e))) --e;
        return e + shift;
    }
}

/// 2^k in T, for k inside T's exponent range, even beyond double's.
template <typename T>
T pow2(int k) {
    return scale_pow2(T(1), k);
}

} // namespace sw::mp_numerics
