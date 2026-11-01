#pragma once
// mp-numerics -- plane rotations computed four ways, for any number type.
//
// A Givens rotation G = [c s; -s c] takes (a, b) to (r, 0). How r is formed
// decides where it fails:
//
//   givens_unscaled    r = sqrt(a*a + b*b), as MTL5's apply_givens_rotation
//                      (operation/givens.hpp, used by GMRES) formed it up to
//                      5.21.0. a*a overflows once |a| passes sqrt(max) and
//                      underflows below sqrt(min), long before r would.
//   givens_max_scaled  divides by m = max(|a|, |b|) first, as LAPACK's lapy2
//                      and MTL5's householder() up to 5.21.0 do. Costs two
//                      divisions.
//   givens_pow2_scaled scales by an exact power of two so the squares land
//                      near 1; not a LAPACK form, but the remedy for tapered
//                      types that mp-numerics#13 tests. MTL5 adopted it in
//                      5.22.0 as mtl::math::givens_rotation (mtl5#673).
//   givens_safe_scaled E. Anderson, "Algorithm 978: Safe Scaling in the Level 1
//                      BLAS", ACM TOMS 44(1), 2017; Reference-LAPACK's dlartg.f90
//                      since 3.10. Scales only when an entry lies outside
//                      [rtmin, rtmax], thresholds derived from the type's range.
//
// The thresholds of the safe form come from std::numeric_limits<T>, as dlartg
// takes them from the Fortran model numbers. For an IEEE type they keep every
// square in range; whether they mean the same thing for a tapered type (posit,
// takum), whose precision falls off toward the ends of its range, is one of
// the questions mp-numerics#13 measures.

#include <algorithm>
#include <cmath>
#include <limits>

#include <sw/mp_numerics/pow2.hpp>

namespace sw::mp_numerics {

template <typename T>
struct plane_rotation {
    T c, s, r;
};

/// r = sqrt(a*a + b*b), c = a/r, s = b/r, as MTL5's givens.hpp computed it up
/// to 5.21.0.
template <typename T>
plane_rotation<T> givens_unscaled(const T& a, const T& b) {
    using std::sqrt;
    const T r = sqrt(a * a + b * b);
    if (r == T(0)) return {T(1), T(0), T(0)};
    return {a / r, b / r, r};
}

/// lapy2's scaling: r = m * sqrt((a/m)^2 + (b/m)^2), m = max(|a|, |b|).
template <typename T>
plane_rotation<T> givens_max_scaled(const T& a, const T& b) {
    using std::abs;
    using std::sqrt;
    const T aa = abs(a), ab = abs(b);
    const T m = aa > ab ? aa : ab;
    if (m == T(0)) return {T(1), T(0), T(0)};
    const T a1 = a / m, b1 = b / m;
    const T r = m * sqrt(a1 * a1 + b1 * b1);
    return {a / r, b / r, r};
}

/// Anderson's safe scaling, transcribed from Reference-LAPACK dlartg.f90
/// (3.10). r carries the sign of a, so c >= 0; G (a, b) = (r, 0) still holds.
template <typename T>
plane_rotation<T> givens_safe_scaled(const T& f, const T& g) {
    using std::abs;
    using std::sqrt;
    using L = std::numeric_limits<T>;
    // safmin = radix^max(minexponent - 1, 1 - maxexponent), from the model
    // numbers as dlartg takes them: the smallest normal, raised if need be so
    // that 1/safmin = safmax is itself a value of T. Forming it as
    // max(min(), 1/max()) instead is not the same: in a fixed-point type
    // 1/max() rounds down, and 1/safmin then overflows.
    static_assert(L::radix == 2, "givens_safe_scaled assumes a binary type");
    const T safmin = pow2<T>(std::max(L::min_exponent - 1, 1 - L::max_exponent));
    const T safmax = pow2<T>(-std::max(L::min_exponent - 1, 1 - L::max_exponent));
    const T rtmin = sqrt(safmin);
    const T rtmax = sqrt(safmax / T(2));

    const T f1 = abs(f), g1 = abs(g);
    if (g == T(0)) return {T(1), T(0), f};
    if (f == T(0)) return {T(0), g > T(0) ? T(1) : T(-1), g1};
    if (f1 > rtmin && f1 < rtmax && g1 > rtmin && g1 < rtmax) {
        const T d = sqrt(f * f + g * g);
        const T r = f < T(0) ? T(-d) : d;
        return {f1 / d, g / r, r};
    }
    T u = f1 > g1 ? f1 : g1;
    if (u < safmin) u = safmin;
    if (u > safmax) u = safmax;
    const T fs = f / u, gs = g / u;
    const T d = sqrt(fs * fs + gs * gs);
    const T r = f < T(0) ? T(-d) : d;
    return {abs(fs) / d, gs / r, r * u};
}

/// Scaled by an exact power of two: with 2^e <= max(|a|, |b|) < 2^(e+1), the
/// scaled entries are at most 2 in magnitude and the largest is at least 1,
/// so the squares sit near 1, where a tapered type (posit, takum) has its
/// most bits, and the scaling itself rounds nothing in a floating-point type.
/// c and s come from the scaled values, as in givens_safe_scaled. Not a
/// LAPACK form: the remedy mp-numerics#13 tests for tapered types.
template <typename T>
plane_rotation<T> givens_pow2_scaled(const T& a, const T& b) {
    using std::abs;
    using std::sqrt;
    const T aa = abs(a), ab = abs(b);
    const T m = aa > ab ? aa : ab;
    if (m == T(0)) return {T(1), T(0), T(0)};
    const int e = binary_exponent(m);
    const T a1 = scale_pow2(a, -e), b1 = scale_pow2(b, -e);
    const T d = sqrt(a1 * a1 + b1 * b1);
    return {a1 / d, b1 / d, scale_pow2(d, e)};
}

/// (u, w) -> (c u + s w, c w - s u), rounded as givens.hpp applies it.
template <typename T>
void apply_rotation(const plane_rotation<T>& g, T& u, T& w) {
    const T t = g.c * u + g.s * w;
    w = g.c * w - g.s * u;
    u = t;
}

} // namespace sw::mp_numerics
