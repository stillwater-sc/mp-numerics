// test_householder_qr: qr_factor_stored_beta (householder_qr.hpp), the
// Householder QR that stores beta on the diagonal as dgeqr2 does, against
// MTL5's qr_factor. Returns non-zero on failure.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>

#include <mtl/mat/dense2D.hpp>
#include <mtl/operation/householder.hpp>
#include <mtl/operation/qr.hpp>
#include <mtl/vec/dense_vector.hpp>

#include <sw/mp_numerics/householder_qr.hpp>

#include <universal/number/fixpnt/fixpnt.hpp>
#include <universal/number/posit/posit.hpp>

namespace {

using namespace sw::mp_numerics;

int failures = 0;

void check(bool ok, const char* what) {
    if (!ok) { std::cerr << "FAIL: " << what << '\n'; ++failures; }
    else std::cout << "ok:   " << what << '\n';
}

template <typename T>
mtl::vec::dense_vector<T> random_vector(std::size_t n, std::mt19937_64& g) {
    mtl::vec::dense_vector<T> x(n);
    for (std::size_t i = 0; i < n; ++i) x(i) = T(double(g() >> 11) * 0x1.0p-53 - 0.5);
    return x;
}

/// householder_beta is the beta mtl::householder forms: from it, dlarfg's
/// v0 = x0 - beta (scaled), tau = -v0 / beta and v(i) = x(i) / v0 rebuild
/// MTL5's tau and v bit for bit. All in the scaled variables, as there.
template <typename T>
bool same_reflector(int trials) {
    using std::abs;
    using scaling = mtl::math::pow2_scaling<T>;
    std::mt19937_64 g(29);
    for (int t = 0; t < trials; ++t) {
        const auto x = random_vector<T>(2 + std::size_t(t % 7), g);
        const auto [v, tau] = mtl::householder(x);
        if (tau == T(0)) continue;
        T m(0);
        for (std::size_t i = 0; i < x.size(); ++i)
            if (abs(x(i)) > m) m = abs(x(i));
        const int e = scaling::exponent(m);
        const T beta = scaling::scale(householder_beta(x), -e);
        const T v0 = scaling::scale(x(0), -e) - beta;
        if (!(tau == T(-v0 / beta))) return false;
        for (std::size_t i = 1; i < x.size(); ++i)
            if (!(v(i) == T(scaling::scale(x(i), -e) / v0))) return false;
    }
    return true;
}

} // namespace

int main() {
    check(same_reflector<double>(1000), "double: householder_beta rebuilds mtl::householder's v and tau");
    check(same_reflector<sw::universal::posit<32, 2>>(1000),
          "posit<32,2>: householder_beta rebuilds mtl::householder's v and tau");

    // A fixed-point column at resolution: x = (0, u, u, u, u), u = 2^-8, so
    // ||x|| = 2u is representable. The products v(i) x(i) = u/2 round to 0,
    // so recomputing R(0,0) as x0 - tau (v^T x) gives 0; beta is -2u.
    {
        using F = sw::universal::fixpnt<16, 8, sw::universal::Saturate, std::uint16_t>;
        const F u(std::ldexp(1.0, -8));
        mtl::mat::dense2D<F> A(5, 1);
        A(0, 0) = F(0);
        for (std::size_t i = 1; i < 5; ++i) A(i, 0) = u;
        mtl::mat::dense2D<F> B(A);
        mtl::vec::dense_vector<F> tau(1), tb(1);
        mtl::qr_factor(A, tau);
        qr_factor_stored_beta(B, tb);
        std::cout << "      fixpnt<16,8> (0, u, u, u, u): qr_factor R(0,0) = " << double(A(0, 0))
                  << ", stored beta = " << double(B(0, 0)) << '\n';
        check(B(0, 0) == F(-2.0 * double(u)), "fixpnt<16,8>: the pivot is beta = -2u, not 0");
    }

    // In double the two differ only by the rounding of the diagonal: the rest
    // of R, and Q, follow from it.
    {
        std::mt19937_64 g(31);
        const std::size_t m = 12, n = 8;
        mtl::mat::dense2D<double> A(m, n);
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < n; ++j) A(i, j) = double(g() >> 11) * 0x1.0p-53 - 0.5;
        mtl::mat::dense2D<double> F(A), G(A);
        mtl::vec::dense_vector<double> tf(n), tg(n);
        mtl::qr_factor(F, tf);
        qr_factor_stored_beta(G, tg);
        const auto Rf = mtl::qr_extract_R(F), Rg = mtl::qr_extract_R(G);
        const auto Qg = mtl::qr_extract_Q(G, tg);
        double dr = 0, res = 0;
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                dr = std::max(dr, std::abs(Rf(i, j) - Rg(i, j)));
                double s = A(i, j);
                for (std::size_t k = 0; k < m; ++k) s -= Qg(i, k) * Rg(k, j);
                res = std::max(res, std::abs(s));
            }
        const double eps = std::numeric_limits<double>::epsilon();
        std::cout << "      double 12x8: max |R - R_stored| = " << dr / eps << " eps, max |A - QR| = " << res / eps
                  << " eps\n";
        check(dr <= 16 * eps, "double: R agrees with qr_factor's to rounding");
        check(res <= 16 * eps, "double: the packed result reads back through qr_extract_Q");
    }

    return failures == 0 ? 0 : 1;
}
