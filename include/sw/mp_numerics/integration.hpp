#pragma once
// mp-numerics -- numerical quadrature (mixed precision).
//
// Composite Newton-Cotes rules, generic over the arithmetic type so the same
// kernel runs in double, float, posit, cfloat, ... and the study can observe how
// the working precision limits the achievable quadrature accuracy. These are
// the composition-layer home for the mixed-precision integration experiments;
// they use only scalar arithmetic on the number type.

#include <cstddef>

namespace sw::mp_numerics {

/// Composite trapezoidal rule for the integral of f over [a, b] with n
/// subintervals. O(h^2) accurate.
template <typename Real, typename F>
Real composite_trapezoidal(F f, Real a, Real b, std::size_t n) {
    if (n < 1) n = 1;
    const Real h = (b - a) / Real(static_cast<long long>(n));
    Real s = (f(a) + f(b)) / Real(2);
    for (std::size_t i = 1; i < n; ++i) s = s + f(a + Real(static_cast<long long>(i)) * h);
    return h * s;
}

/// Composite Simpson's rule for the integral of f over [a, b] with n
/// subintervals (n is rounded up to the next even value). O(h^4) accurate, and
/// exact for cubics.
template <typename Real, typename F>
Real composite_simpson(F f, Real a, Real b, std::size_t n) {
    if (n < 2) n = 2;
    if (n % 2 != 0) ++n;                                   // Simpson needs an even count
    const Real h = (b - a) / Real(static_cast<long long>(n));
    Real s = f(a) + f(b);
    for (std::size_t i = 1; i < n; ++i) {
        const Real x = a + Real(static_cast<long long>(i)) * h;
        s = s + (i % 2 != 0 ? Real(4) : Real(2)) * f(x);
    }
    return (h / Real(3)) * s;
}

} // namespace sw::mp_numerics
