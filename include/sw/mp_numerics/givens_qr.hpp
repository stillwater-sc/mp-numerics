#pragma once
// mp-numerics -- QR by Givens rotations, for any number type.
//
// MTL5 factors QR by Householder reflectors only (operation/qr.hpp). This is
// the Givens form the decompositions study (mp-numerics#14) compares them
// with. Each rotation comes from MTL5's own mtl::math::givens_rotation
// (math/pow2_scaling.hpp, 5.22.0), which scales a and b by an exact power of
// two before squaring, so a difference from qr_factor is the method's and not
// the rotation's.
//
// Column by column, bottom up (Golub & Van Loan, 4th ed., Alg. 5.2.4): entry
// (i, j) is zeroed against (i-1, j), and the rotation is applied to rows i-1
// and i of the trailing columns and of whatever else is carried along (Q^T
// for givens_qr, b for givens_qr_solve). An entry that is already zero is
// skipped, so a banded or Hessenberg A keeps its zeros.

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <utility>

#include <mtl/mat/dense2D.hpp>
#include <mtl/math/pow2_scaling.hpp>
#include <mtl/vec/dense_vector.hpp>

namespace sw::mp_numerics {

/// A = Q R with Q (m x m) orthogonal and R (m x n) upper triangular.
template <typename T>
struct givens_qr_result {
    mtl::mat::dense2D<T> Q, R;
};

namespace detail {

/// Reduce R (m x n, m >= n) to upper triangular form in place. For every
/// rotation, carry(i, c, s) applies [c s; -s c] to rows i-1 and i of what
/// the caller carries along.
template <typename T, typename Carry>
void givens_triangularize(mtl::mat::dense2D<T>& R, Carry&& carry) {
    const std::size_t m = R.num_rows(), n = R.num_cols();
    assert(m >= n);
    if (m == 0) return;
    const T zero(0);
    for (std::size_t j = 0; j < std::min(n, m - 1); ++j) {
        for (std::size_t i = m - 1; i > j; --i) {
            if (R(i, j) == zero) continue;
            const auto g = mtl::math::givens_rotation(R(i - 1, j), R(i, j));
            R(i - 1, j) = g.r;
            R(i, j) = zero;
            for (std::size_t col = j + 1; col < n; ++col) {
                const T x = R(i - 1, col), y = R(i, col);
                R(i - 1, col) = g.c * x + g.s * y;
                R(i, col) = g.c * y - g.s * x;
            }
            carry(i, g.c, g.s);
        }
    }
}

} // namespace detail

/// Q and R of A (m x n, m >= n) by Givens rotations.
template <typename T>
givens_qr_result<T> givens_qr(const mtl::mat::dense2D<T>& A) {
    const std::size_t m = A.num_rows();
    mtl::mat::dense2D<T> R(A);
    mtl::mat::dense2D<T> Qt(m, m);   // the product of the rotations, G_k ... G_1
    for (std::size_t r = 0; r < m; ++r)
        for (std::size_t c = 0; c < m; ++c) Qt(r, c) = r == c ? T(1) : T(0);
    detail::givens_triangularize(R, [&](std::size_t i, const T& c, const T& s) {
        for (std::size_t col = 0; col < m; ++col) {
            const T x = Qt(i - 1, col), y = Qt(i, col);
            Qt(i - 1, col) = c * x + s * y;
            Qt(i, col) = c * y - s * x;
        }
    });
    mtl::mat::dense2D<T> Q(m, m);
    for (std::size_t r = 0; r < m; ++r)
        for (std::size_t c = 0; c < m; ++c) Q(r, c) = Qt(c, r);
    return {std::move(Q), std::move(R)};
}

/// The least-squares solution of A x = b (A m x n, m >= n, full column rank):
/// Q^T b by the same rotations that triangularize A, then back substitution
/// with R's leading n x n block.
template <typename T>
mtl::vec::dense_vector<T> givens_qr_solve(const mtl::mat::dense2D<T>& A, const mtl::vec::dense_vector<T>& b) {
    const std::size_t m = A.num_rows(), n = A.num_cols();
    assert(b.size() == m);
    mtl::mat::dense2D<T> R(A);
    mtl::vec::dense_vector<T> y(m);
    for (std::size_t i = 0; i < m; ++i) y(i) = b(i);
    detail::givens_triangularize(R, [&](std::size_t i, const T& c, const T& s) {
        const T u = y(i - 1), w = y(i);
        y(i - 1) = c * u + s * w;
        y(i) = c * w - s * u;
    });
    mtl::vec::dense_vector<T> x(n);
    for (std::size_t ii = 0; ii < n; ++ii) {
        const std::size_t i = n - 1 - ii;
        T sum(0);
        for (std::size_t j = i + 1; j < n; ++j) sum += R(i, j) * x(j);
        x(i) = (y(i) - sum) / R(i, i);
    }
    return x;
}

} // namespace sw::mp_numerics
