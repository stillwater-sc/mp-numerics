#pragma once
// mp-numerics -- scalar root finding (mixed precision).
//
// Bisection, Newton, and secant iterations, generic over the arithmetic type so
// the same solver runs in double, float, posit, cfloat, ... and the study can
// observe how the working precision limits the achievable root accuracy. Scalar
// arithmetic only; the heavier polynomial root finders (Jenkins-Traub rpoly) are
// a follow-up (see docs/roadmap.md).

#include <cstddef>

namespace sw::mp_numerics {

/// Outcome of a root-finding iteration.
struct root_result {
    std::size_t iterations = 0;
    bool        converged  = false;   ///< reached |f| or |step| <= tol
};

namespace detail {
    template <typename T> T mag(const T& x) { return x < T(0) ? T(-x) : x; }
}

/// Bisection on [a, b] (requires a sign change f(a)*f(b) <= 0). Robust, linear.
template <typename Real, typename F>
Real bisection(F f, Real a, Real b, Real tol, std::size_t maxit, root_result* out = nullptr) {
    Real fa = f(a), m = a;
    root_result r;
    for (r.iterations = 0; r.iterations < maxit; ++r.iterations) {
        m = (a + b) / Real(2);
        Real fm = f(m);
        if (detail::mag(fm) <= tol || (detail::mag(b - a) / Real(2)) <= tol) { r.converged = true; break; }
        if ((fa < Real(0)) == (fm < Real(0))) { a = m; fa = fm; } else { b = m; }
    }
    if (out) *out = r;
    return m;
}

/// Newton's method from x0 using f and its derivative df. Quadratic near a
/// simple root; falls back to no progress if df underflows to zero.
template <typename Real, typename F, typename DF>
Real newton(F f, DF df, Real x0, Real tol, std::size_t maxit, root_result* out = nullptr) {
    Real x = x0;
    root_result r;
    for (r.iterations = 0; r.iterations < maxit; ++r.iterations) {
        Real fx = f(x);
        if (detail::mag(fx) <= tol) { r.converged = true; break; }
        Real d = df(x);
        if (detail::mag(d) == Real(0)) break;         // derivative underflow -> stop
        Real step = fx / d;
        x = x - step;
        if (detail::mag(step) <= tol) { r.converged = detail::mag(f(x)) <= tol; break; }
    }
    if (out) *out = r;
    return x;
}

/// Secant method from x0, x1 (no derivative needed). Superlinear.
template <typename Real, typename F>
Real secant(F f, Real x0, Real x1, Real tol, std::size_t maxit, root_result* out = nullptr) {
    Real a = x0, b = x1, fa = f(a), fb = f(b);
    root_result r;
    for (r.iterations = 0; r.iterations < maxit; ++r.iterations) {
        if (detail::mag(fb) <= tol) { r.converged = true; break; }
        Real denom = fb - fa;
        if (detail::mag(denom) == Real(0)) break;
        Real x = b - fb * (b - a) / denom;
        a = b; fa = fb; b = x; fb = f(b);
        if (detail::mag(b - a) <= tol) { r.converged = detail::mag(fb) <= tol; break; }
    }
    if (out) *out = r;
    return b;
}

} // namespace sw::mp_numerics
