#pragma once
// mp-numerics -- one-dimensional optimization (mixed precision).
//
// Derivative-free golden-section minimization, generic over the arithmetic type
// so the same search runs in double, float, posit, cfloat, ... and the study can
// observe how the working precision limits the locatable minimizer accuracy
// (the bracket cannot shrink below the type's spacing near the minimum). Scalar
// arithmetic only.

#include <cstddef>

namespace sw::mp_numerics {

/// Outcome of a 1-D minimization.
struct min_result {
    std::size_t iterations = 0;
    bool        converged  = false;   ///< bracket width <= tol
};

namespace detail {
    template <typename T> T amag(const T& x) { return x < T(0) ? T(-x) : x; }
}

/// Golden-section search for a minimizer of a unimodal f on [a, b]. Returns the
/// bracket midpoint; sets `out` if provided.
template <typename Real, typename F>
Real golden_section_min(F f, Real a, Real b, Real tol, std::size_t maxit, min_result* out = nullptr) {
    const Real gr = Real(0.6180339887498949);          // (sqrt(5) - 1) / 2
    Real c = b - gr * (b - a);
    Real d = a + gr * (b - a);
    Real fc = f(c), fd = f(d);
    min_result r;
    for (r.iterations = 0; r.iterations < maxit; ++r.iterations) {
        if (detail::amag(b - a) <= tol) { r.converged = true; break; }
        if (fc < fd) { b = d; d = c; fd = fc; c = b - gr * (b - a); fc = f(c); }
        else         { a = c; c = d; fc = fd; d = a + gr * (b - a); fd = f(d); }
    }
    if (out) *out = r;
    return (a + b) / Real(2);
}

} // namespace sw::mp_numerics
