// test_givens_qr: the Givens QR reference (givens_qr.hpp) on cases whose
// outcome is known, against MTL5's Householder qr_factor. Returns non-zero on
// failure.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>
#include <utility>

#include <mtl/mat/dense2D.hpp>
#include <mtl/operation/qr.hpp>
#include <mtl/vec/dense_vector.hpp>

#include <sw/mp_numerics/givens_qr.hpp>

#include <universal/number/posit/posit.hpp>

namespace {

using namespace sw::mp_numerics;

int failures = 0;

void check(bool ok, const char* what) {
    if (!ok) { std::cerr << "FAIL: " << what << '\n'; ++failures; }
    else std::cout << "ok:   " << what << '\n';
}

template <typename T>
mtl::mat::dense2D<T> random_matrix(std::size_t m, std::size_t n, std::uint64_t seed, double scale = 1.0) {
    std::mt19937_64 g(seed);
    mtl::mat::dense2D<T> A(m, n);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < n; ++j) A(i, j) = T(scale * (double(g() >> 11) * 0x1.0p-53 - 0.5));
    return A;
}

/// ||A - Q R||_F / ||A||_F, in double.
template <typename T>
double residual(const mtl::mat::dense2D<T>& A, const mtl::mat::dense2D<T>& Q, const mtl::mat::dense2D<T>& R) {
    double num = 0, den = 0;
    for (std::size_t i = 0; i < A.num_rows(); ++i)
        for (std::size_t j = 0; j < A.num_cols(); ++j) {
            double s = double(A(i, j));
            for (std::size_t k = 0; k < Q.num_cols(); ++k) s -= double(Q(i, k)) * double(R(k, j));
            num += s * s;
            den += double(A(i, j)) * double(A(i, j));
        }
    return std::sqrt(num / den);
}

/// ||Q^T Q - I||_F, in double.
template <typename T>
double orthogonality(const mtl::mat::dense2D<T>& Q) {
    double sum = 0;
    for (std::size_t i = 0; i < Q.num_cols(); ++i)
        for (std::size_t j = 0; j < Q.num_cols(); ++j) {
            double s = i == j ? -1.0 : 0.0;
            for (std::size_t k = 0; k < Q.num_rows(); ++k) s += double(Q(k, i)) * double(Q(k, j));
            sum += s * s;
        }
    return std::sqrt(sum);
}

template <typename T>
bool upper_triangular(const mtl::mat::dense2D<T>& R) {
    for (std::size_t i = 0; i < R.num_rows(); ++i)
        for (std::size_t j = 0; j < std::min(i, R.num_cols()); ++j)
            if (R(i, j) != T(0)) return false;
    return true;
}

} // namespace

int main() {
    const double eps = std::numeric_limits<double>::epsilon();

    // A square and a tall double matrix: backward stable, orthogonal Q,
    // exactly triangular R.
    for (auto [m, n] : {std::pair<std::size_t, std::size_t>{8, 8}, {16, 8}}) {
        const auto A = random_matrix<double>(m, n, 7 + m);
        const auto [Q, R] = givens_qr(A);
        std::cout << "      double " << m << "x" << n << ": residual " << residual(A, Q, R) / eps
                  << " eps, orthogonality " << orthogonality(Q) / eps << " eps\n";
        check(residual(A, Q, R) <= 8 * double(m) * eps, "double: ||A - QR|| / ||A|| within 8 m eps");
        check(orthogonality(Q) <= 8 * double(m) * eps, "double: ||Q^T Q - I|| within 8 m eps");
        check(upper_triangular(R), "double: R is exactly upper triangular");
    }

    // With R's diagonal made positive, the QR of a full-rank A is unique, so
    // Givens and Householder must agree to rounding.
    {
        const auto A = random_matrix<double>(8, 8, 11);
        auto [Qg, Rg] = givens_qr(A);
        mtl::mat::dense2D<double> F(A);
        mtl::vec::dense_vector<double> tau(8);
        mtl::qr_factor(F, tau);
        auto Qh = mtl::qr_extract_Q(F, tau);
        auto Rh = mtl::qr_extract_R(F);
        mtl::qr_normalize_signs(Qg, Rg);
        mtl::qr_normalize_signs(Qh, Rh);
        double d = 0;
        for (std::size_t i = 0; i < 8; ++i)
            for (std::size_t j = 0; j < 8; ++j) d = std::max(d, std::abs(Rg(i, j) - Rh(i, j)));
        std::cout << "      max |R_givens - R_householder| = " << d / eps << " eps\n";
        check(d <= 64 * eps, "Givens and Householder give the same R up to signs");
    }

    // Least squares: a consistent tall system recovers x.
    {
        const auto A = random_matrix<double>(16, 8, 13);
        mtl::vec::dense_vector<double> x0(8), b(16);
        for (std::size_t j = 0; j < 8; ++j) x0(j) = double(j + 1) / 8;
        for (std::size_t i = 0; i < 16; ++i) {
            b(i) = 0;
            for (std::size_t j = 0; j < 8; ++j) b(i) += A(i, j) * x0(j);
        }
        const auto x = givens_qr_solve(A, b);
        double d = 0;
        for (std::size_t j = 0; j < 8; ++j) d = std::max(d, std::abs(x(j) - x0(j)));
        check(d <= 1e-12, "givens_qr_solve recovers x of a consistent system");
    }

    // A zero entry is skipped, so a Hessenberg A takes one rotation per column
    // and its zeros below the subdiagonal are never touched.
    {
        auto A = random_matrix<double>(6, 6, 17);
        for (std::size_t i = 0; i < 6; ++i)
            for (std::size_t j = 0; j + 1 < i; ++j) A(i, j) = 0;
        const auto [Q, R] = givens_qr(A);
        bool hessenberg_q = true;   // Q of a Hessenberg A is Hessenberg too
        for (std::size_t i = 0; i < 6; ++i)
            for (std::size_t j = 0; j + 1 < i; ++j) hessenberg_q = hessenberg_q && Q(i, j) == 0.0;
        check(hessenberg_q && residual(A, Q, R) <= 64 * eps, "Hessenberg A: Q stays Hessenberg");
    }

    // The rotations are scaled by an exact power of two (MTL5 5.22.0), so an
    // A with entries near 2^600, whose squares overflow double, still factors.
    {
        const auto A = random_matrix<double>(8, 8, 19, std::ldexp(1.0, 600));
        const auto [Q, R] = givens_qr(A);
        bool finite = true;
        for (std::size_t i = 0; i < 8; ++i)
            for (std::size_t j = 0; j < 8; ++j) finite = finite && std::isfinite(Q(i, j)) && std::isfinite(R(i, j));
        check(finite && orthogonality(Q) <= 64 * eps, "double near 2^600: finite and orthogonal");
    }

    // It is generic in the element type.
    {
        using P = sw::universal::posit<32, 2>;
        const auto A = random_matrix<P>(8, 8, 23);
        const auto [Q, R] = givens_qr(A);
        const double ep = double(std::numeric_limits<P>::epsilon());
        std::cout << "      posit<32,2>: residual " << residual(A, Q, R) / ep << " eps, orthogonality "
                  << orthogonality(Q) / ep << " eps\n";
        check(residual(A, Q, R) <= 64 * ep && orthogonality(Q) <= 64 * ep, "posit<32,2>: backward stable");
    }

    return failures == 0 ? 0 : 1;
}
