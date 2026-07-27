#pragma once
// mp-numerics -- real-coefficient polynomial root finding (Jenkins-Traub, "rpoly").
//
// All roots of a degree-n polynomial with real coefficients via the three-stage
// Jenkins-Traub algorithm (M. Jenkins & J. Traub, "A Three-Stage Algorithm for
// Real Polynomials Using Quadratic Iteration", SIAM J. Numer. Anal. 1970).
// Generic over the arithmetic type `Real`, so the same solver runs in double,
// float, posit, cfloat, ... and the precision study can observe how the working
// precision limits achievable root accuracy.
//
// The engine is the canonical ak1 / Netlib TOMS-493 reference implementation,
// templatized over `Real` in detail/rpoly_ak1.hpp (mp-numerics#7). It replaces
// the earlier lower-fidelity port, which lost ~13 orders of accuracy on
// well-separated roots (double: 1.6e-2 vs the ak1 engine's 3.6e-15). One key
// adaptation makes it correct for tapered-precision types: the classic algorithm
// unconditionally power-of-2 rescales the coefficients, which pushes O(1)
// coefficients into a posit's low-precision tapered region; detail/rpoly_ak1.hpp
// only rescales when the coefficients actually approach the type's representable
// extremes, so every number type now reaches its precision floor (posit<64,2> is
// exact on well-separated integer roots).
//
// Coefficients are highest-degree first: {a_n, ..., a_1, a_0} for
// a_n*x^n + ... + a_1*x + a_0. Complex roots use sw::universal::complex<Real>
// (std::complex is undefined behavior for the Universal number types).

#include <cstddef>
#include <vector>

#include <mtl/vec/dense_vector.hpp>
#include <universal/math/complex.hpp>

#include <sw/mp_numerics/detail/rpoly_ak1.hpp>

namespace sw::mp_numerics {

template <typename Real>
using Complex = sw::universal::complex<Real>;

/// All n roots of a real polynomial, plus a success flag.
template <typename Real>
struct poly_roots_result {
    std::vector<Complex<Real>> roots;   ///< the roots (complex conjugate pairs for complex roots)
    bool ok = false;                    ///< false if the leading coefficient was zero / extraction failed / degree out of range
};

/// All roots of a real-coefficient polynomial (highest-degree-first coefficients)
/// via Jenkins-Traub, returned as complex values in the working precision `Real`.
/// Supports degree 1..detail::ak1::MAXDEGREE; out-of-range input returns ok=false.
template <typename Real>
[[nodiscard]] poly_roots_result<Real> poly_roots(const mtl::vec::dense_vector<Real>& coefficients) {
    poly_roots_result<Real> out;
    const int degree = static_cast<int>(coefficients.size()) - 1;
    if (degree < 1 || degree > detail::ak1::MAXDEGREE) return out;   // ok stays false

    Real op[detail::ak1::MDP1]{};
    for (int i = 0; i <= degree; ++i) op[i] = coefficients[i];

    Real zr[detail::ak1::MAXDEGREE]{}, zi[detail::ak1::MAXDEGREE]{};
    int n = degree;
    detail::ak1::rpoly_ak1<Real>(op, &n, zr, zi);

    // ak1 sets *Degree to the number of roots recovered (degree on success, a
    // smaller count on partial convergence, 0/-1 on failure).
    if (n < 1) return out;                                           // ok stays false
    out.ok = true;
    out.roots.reserve(n);
    for (int i = 0; i < n; ++i) out.roots.emplace_back(zr[i], zi[i]);
    return out;
}

/// Build a monic polynomial (highest-degree-first) from its real roots: prod (x - r_i).
template <typename Real>
[[nodiscard]] mtl::vec::dense_vector<Real> poly_from_real_roots(const std::vector<Real>& roots) {
    std::vector<Real> p{ Real(1) };                                 // p(x) = 1
    for (const Real& r : roots) {                                   // multiply by (x - r)
        std::vector<Real> q(p.size() + 1, Real(0));
        for (std::size_t i = 0; i < p.size(); ++i) {
            q[i]     += p[i];                                        // x * p(x)
            q[i + 1] += -r * p[i];                                  // -r * p(x)
        }
        p.swap(q);
    }
    mtl::vec::dense_vector<Real> out(p.size());
    for (std::size_t i = 0; i < p.size(); ++i) out[i] = p[i];
    return out;
}

} // namespace sw::mp_numerics
