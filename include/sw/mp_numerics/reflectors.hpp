#pragma once
// mp-numerics -- a Householder reflector in LAPACK dlarfg's form, scaled by an
// exact power of two, for any number type.
//
// MTL5's householder() (operation/householder.hpp) reflects x onto +||x|| e1
// when x0 > 0 (Parlett's choice of sign) and normalizes v so that v(0) = 1.
// With that sign v0 = -sigma / (x0 + ||x||) is small when the tail of x is
// small, and v1 = x1 / v0 ~ -2 x0 / x1 is large: in a tapered type (posit,
// takum) that is where the type has the fewest bits. dlarfg reflects onto
// beta e1 with beta = -sign(x0) ||x||, so |v0| = |x0| + ||x|| >= ||x|| and
// every entry of the normalized v is at most 1 in magnitude.
//
// householder_lapack_pow2 is that form, with x first scaled by 2^-e,
// 2^e <= max |x_i|. The scaling is exact in a floating-point type and puts
// the sum of squares near 1, inside every type's range and where a tapered
// type is most precise. It is not a transcription of dlarfg, which scales
// only when the norm under- or overflows and forms v by multiplying by a
// reciprocal; mp-numerics#13 tests it as the reflector for tapered types.

#include <cmath>

#include <mtl/vec/dense_vector.hpp>

#include <sw/mp_numerics/pow2.hpp>

namespace sw::mp_numerics {

/// H = I - tau v v^T with v(0) = 1 and H x = beta e1. tau == 0 is the identity.
template <typename T>
struct normalized_reflector {
    mtl::vec::dense_vector<T> v;
    T tau;
    T beta;
};

template <typename T>
normalized_reflector<T> householder_lapack_pow2(const mtl::vec::dense_vector<T>& x) {
    using std::abs;
    using std::sqrt;
    const std::size_t n = x.size();
    normalized_reflector<T> h{mtl::vec::dense_vector<T>(n), T(0), n ? T(x(0)) : T(0)};
    for (std::size_t i = 0; i < n; ++i) h.v(i) = T(0);
    if (n == 0) return h;
    h.v(0) = T(1);

    T m(0);
    for (std::size_t i = 1; i < n; ++i)
        if (abs(x(i)) > m) m = abs(x(i));
    if (m == T(0)) return h;   // the tail is zero: H = I, beta = x0
    if (abs(x(0)) > m) m = abs(x(0));
    const int e = binary_exponent(m);

    mtl::vec::dense_vector<T> xs(n);
    for (std::size_t i = 0; i < n; ++i) xs(i) = scale_pow2(T(x(i)), -e);
    T sum(0);
    for (std::size_t i = 0; i < n; ++i) sum += xs(i) * xs(i);
    const T norm = sqrt(sum);
    const T alpha = xs(0);
    const T beta = alpha >= T(0) ? T(-norm) : norm;
    const T d = alpha - beta;   // |d| = |alpha| + norm, no cancellation
    h.tau = (beta - alpha) / beta;
    for (std::size_t i = 1; i < n; ++i) h.v(i) = xs(i) / d;
    h.beta = scale_pow2(beta, e);
    return h;
}

} // namespace sw::mp_numerics
