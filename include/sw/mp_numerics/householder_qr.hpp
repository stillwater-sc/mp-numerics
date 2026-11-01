#pragma once
// mp-numerics -- Householder QR that stores each reflector's beta on the
// diagonal, as LAPACK's dgeqr2 does.
//
// MTL5's generic qr_factor (operation/qr.hpp, 5.22.0) applies column j's
// reflector to columns j onward, column j included, so R(j,j) is recomputed
// as x0 - tau (v^T x) rather than set to beta, where H x = beta e1. In exact
// arithmetic the two agree. In fixed point the products v(i) x(i) of a column
// near the type's resolution round to zero, and R(j,j) comes out 0 although
// beta = -sign(x0) ||x|| is representable; qr_solve then divides by zero.
// dgeqr2, and MTL5's own textbook path, store beta and apply the reflector
// to columns j+1 onward only. The decompositions study (mp-numerics#14)
// tests this form as the remedy.
//
// The reflector is MTL5's: mtl::householder gives v and tau, and
// householder_beta forms beta as mtl::householder forms it internally (it
// returns only v and tau). The result is in qr_factor's packed layout, so
// mtl::qr_extract_Q, mtl::qr_extract_R and mtl::qr_solve read it unchanged.

#include <algorithm>
#include <cmath>
#include <cstddef>

#include <mtl/mat/dense2D.hpp>
#include <mtl/math/pow2_scaling.hpp>
#include <mtl/operation/householder.hpp>
#include <mtl/vec/dense_vector.hpp>

namespace sw::mp_numerics {

/// beta of a real x's reflector, H x = beta e1, as mtl::householder (MTL5
/// 5.22.0) forms it: x scaled by an exact 2^-e with 2^e <= max |x_i|, the sum
/// of squares of the tail, beta = -sign(x0) sqrt(x0^2 + sigma) with
/// sign(0) = +1, scaled back. Only meaningful where mtl::householder returns
/// a nonzero tau.
template <typename T>
T householder_beta(const mtl::vec::dense_vector<T>& x) {
    using std::abs;
    using std::sqrt;
    using scaling = mtl::math::pow2_scaling<T>;
    const std::size_t n = x.size();
    T m(0);
    for (std::size_t i = 0; i < n; ++i)
        if (abs(x(i)) > m) m = abs(x(i));
    const int e = scaling::exponent(m);
    T sigma(0);
    for (std::size_t i = 1; i < n; ++i) {
        const T xs = scaling::scale(x(i), -e);
        sigma += xs * xs;
    }
    const T x0 = scaling::scale(x(0), -e);
    const T norm = sqrt(x0 * x0 + sigma);
    return scaling::scale(x0 >= T(0) ? T(-norm) : norm, e);
}

/// A (m x n) is overwritten with R on and above the diagonal and the
/// reflectors below it; tau(k) is reflector k's scalar.
template <typename T>
void qr_factor_stored_beta(mtl::mat::dense2D<T>& A, mtl::vec::dense_vector<T>& tau) {
    const std::size_t m = A.num_rows(), n = A.num_cols(), k = std::min(m, n);
    tau.change_dim(k);
    for (std::size_t j = 0; j < k; ++j) {
        mtl::vec::dense_vector<T> col(m - j);
        for (std::size_t i = 0; i < m - j; ++i) col(i) = A(j + i, j);
        auto [v, t] = mtl::householder(col);
        tau(j) = t;
        if (t == T(0)) continue;   // H = I: the column, and its pivot x0, stay as they are
        mtl::apply_householder_left(A, v, t, j, j + 1);
        A(j, j) = householder_beta(col);
        for (std::size_t i = 1; i < m - j; ++i) A(j + i, j) = v(i);
    }
}

} // namespace sw::mp_numerics
