#pragma once
// mp-numerics -- exact error measures of a computed rotation or reflector.
//
// Every measure here is formed exactly in expansion arithmetic (exact.hpp) from
// the computed values, and rounded once to double. So a measure reports the
// error of the operator the arithmetic produced, not the error of measuring it.
//
// The input x = (a, b) and everything derived from it are read out scaled by
// 2^-k, the caller's choice of binary exponent near ||x||. The operator's own
// coefficients (c, s, tau, v1) are scale-free and read unscaled.
//
// A computed value that, so scaled, falls outside double's range cannot be
// read (to_exact throws unreadable_value). It is off from its true size by a
// factor near 2^1000 or more, which only an operator gone wrong produces (in a
// saturating type, a square clamped at max()). Such an operator is reported
// with every error infinite, so it counts as gross.

#include <cmath>
#include <vector>

#include <sw/mp_numerics/exact.hpp>

namespace sw::mp_numerics {

/// The measures of one operator applied to one x.
struct operator_error {
    bool   finite   = true;   ///< every computed value is finite
    bool   identity = false;  ///< the routine returned the identity
    double orth     = 0;      ///< ||Q^T Q - I||_2 of the computed operator
    double r_err    = 0;      ///< |stored surviving entry - ||x||| / ||x||
    double annih    = 0;      ///< |second entry of Q x, exactly| / ||x||
    double annih_fp = 0;      ///< |second entry of Q x, computed in T| / ||x||
};

namespace detail {

/// The measures of an operator with a value to_exact cannot read: every error infinite.
inline operator_error out_of_reach() {
    operator_error e;
    e.orth = e.r_err = e.annih = e.annih_fp = HUGE_VAL;
    return e;
}

/// (|R| - n) / n for R exact and n^2 = N2 exact, as |R^2 - N2| / (n (|R| + n)):
/// the numerator is exact, so no cancellation is left in the subtraction.
inline double relative_length_error(const exact& R, const exact& N2) {
    exact d = R * R;
    d -= N2;
    const double n = std::sqrt(double(N2));
    return abs_double(d) / (n * (std::abs(double(R)) + n));
}

/// The measures below, reading every value exactly; they throw
/// unreadable_value where a value cannot be read.
template <typename T>
operator_error measure_rotation_unguarded(const T& a, const T& b, int k,
                                          const T& c, const T& s, const T& r,
                                          const T& y0, const T& y1) {
    operator_error e;
    for (double v : {scaled_double(c), scaled_double(s), scaled_double(r, k),
                     scaled_double(y0, k), scaled_double(y1, k)})
        if (!std::isfinite(v)) e.finite = false;
    if (!e.finite) return e;

    const exact A = to_exact(a, k), B = to_exact(b, k);
    const exact C = to_exact(c), S = to_exact(s);
    exact N2 = A * A;
    N2 += B * B;
    const double n = std::sqrt(double(N2));

    // G^T G = (c^2 + s^2) I, so ||G^T G - I||_2 = |c^2 + s^2 - 1|.
    exact o = C * C;
    o += S * S;
    o -= exact(1.0);
    e.orth = abs_double(o);

    e.r_err = detail::relative_length_error(to_exact(r, k), N2);
    exact z = C * B;
    z -= S * A;
    e.annih = abs_double(z) / n;
    e.annih_fp = abs_double(to_exact(y1, k)) / n;
    return e;
}

} // namespace detail

/// A rotation G = [c s; -s c] built for x = (a, b). `r` is the value the
/// algorithm stores in place of a; (y0, y1) is G x computed in T.
template <typename T>
operator_error measure_rotation(const T& a, const T& b, int k,
                                const T& c, const T& s, const T& r,
                                const T& y0, const T& y1) {
    try {
        return detail::measure_rotation_unguarded(a, b, k, c, s, r, y0, y1);
    } catch (const unreadable_value&) {
        return detail::out_of_reach();
    }
}

namespace detail {

template <typename T>
operator_error measure_reflector_unguarded(const T& a, const T& b, int k,
                                           const T& tau, const T& v1,
                                           const T& y0, const T& y1) {
    operator_error e;
    for (double v : {scaled_double(tau), scaled_double(v1),
                     scaled_double(y0, k), scaled_double(y1, k)})
        if (!std::isfinite(v)) e.finite = false;
    if (!e.finite) return e;

    const exact A = to_exact(a, k), B = to_exact(b, k);
    exact N2 = A * A;
    N2 += B * B;
    const double n = std::sqrt(double(N2));

    if (tau == T(0)) {
        e.identity = true;
        e.r_err = detail::relative_length_error(A, N2);
        e.annih = e.annih_fp = abs_double(B) / n;
        return e;
    }

    const exact TAU = to_exact(tau), V1 = to_exact(v1);
    // H^T H - I = (g^2 - 2 g) v v^T / (v^T v) with g = tau v^T v, so its
    // 2-norm is |g| |g - 2|.
    exact vv = V1 * V1;
    vv += exact(1.0);
    const exact g = TAU * vv;
    exact gm2 = g;
    gm2 -= exact(2.0);
    e.orth = abs_double(g) * abs_double(gm2);

    e.r_err = detail::relative_length_error(to_exact(y0, k), N2);
    // (H x)_1 = b - tau v1 (a + v1 b), exactly.
    exact d = V1 * B;
    d += A;
    exact y = TAU * V1 * d;
    exact z = B;
    z -= y;
    e.annih = abs_double(z) / n;
    e.annih_fp = abs_double(to_exact(y1, k)) / n;
    return e;
}

} // namespace detail

/// A reflector H = I - tau v v^T with v = (1, v1), built for x = (a, b).
/// (y0, y1) is H x computed in T; y0 is what a QR factorization stores.
/// tau == 0 is MTL5's identity return.
template <typename T>
operator_error measure_reflector(const T& a, const T& b, int k,
                                 const T& tau, const T& v1,
                                 const T& y0, const T& y1) {
    try {
        return detail::measure_reflector_unguarded(a, b, k, tau, v1, y0, y1);
    } catch (const unreadable_value&) {
        return detail::out_of_reach();
    }
}

namespace detail {

template <typename T>
operator_error measure_textbook_reflector_unguarded(const T& a, const T& b, int k,
                                                    const T& v0, const T& v1, const T& vtv,
                                                    const T& alpha, const T& y1, int kv) {
    operator_error e;
    for (double v : {scaled_double(v0, kv), scaled_double(v1, kv), scaled_double(vtv, 2 * kv),
                     scaled_double(alpha, k), scaled_double(y1, k)})
        if (!std::isfinite(v)) e.finite = false;
    if (!e.finite) return e;

    const exact A = to_exact(a, k), B = to_exact(b, k);
    exact N2 = A * A;
    N2 += B * B;
    const double n = std::sqrt(double(N2));

    if (!(vtv > T(0))) {   // skipped column: the identity
        e.identity = true;
        e.r_err = detail::relative_length_error(A, N2);
        e.annih = e.annih_fp = abs_double(B) / n;
        return e;
    }

    const exact V0 = to_exact(v0, kv), V1 = to_exact(v1, kv), W = to_exact(vtv, 2 * kv);
    // g = 2 v^T v / vtv; ||H^T H - I||_2 = |g| |g - 2| and
    // g - 2 = 2 (v^T v - vtv) / vtv, with the numerator exact.
    exact vv = V0 * V0;
    vv += V1 * V1;
    exact gap = vv;
    gap -= W;
    const double w = double(W);
    const double g = 2 * double(vv) / w;
    e.orth = std::abs(g) * 2 * abs_double(gap) / w;

    e.r_err = detail::relative_length_error(to_exact(alpha, k), N2);
    // (H x)_1 = (b vtv - 2 v1 (v0 a + v1 b)) / vtv, the numerator exact.
    exact d = V0 * A;
    d += V1 * B;
    exact z = B * W;
    z -= exact(2.0) * V1 * d;
    e.annih = abs_double(z) / w / n;
    e.annih_fp = abs_double(to_exact(y1, k)) / n;
    return e;
}

} // namespace detail

/// MTL5's textbook reflector (householder(x, tag::textbook_reflector)):
/// H = I - (2 / vtv) v v^T with v = (v0, v1) unnormalized, applied as
/// y - (2 (v^T y) / vtv) v. vtv is the computed v^T v, so the operator applied
/// is exactly I - (2 / vtv) v v^T for that vtv. `alpha` is what QR stores for
/// R(k,k), and y1 the second entry of H x computed in T. v0, v1 and alpha
/// scale like x (by 2^-k), vtv like its square (by 2^-2k).
/// v may carry its own scale: v0 and v1 are read scaled by 2^-kv and vtv by
/// 2^-2kv. Every measure is homogeneous of degree 0 in (v, vtv), so kv only
/// has to keep those values in double's range; it defaults to k.
template <typename T>
operator_error measure_textbook_reflector(const T& a, const T& b, int k,
                                          const T& v0, const T& v1, const T& vtv,
                                          const T& alpha, const T& y1, int kv) {
    try {
        return detail::measure_textbook_reflector_unguarded(a, b, k, v0, v1, vtv, alpha, y1, kv);
    } catch (const unreadable_value&) {
        return detail::out_of_reach();
    }
}

/// ||Q^T Q - I||_F of an n x n matrix given row-major as values of T, exactly.
template <typename T>
double orthogonality_frobenius(const std::vector<T>& q, int n) {
    std::vector<exact> Q;
    Q.reserve(q.size());
    for (const T& v : q) {
        if (!std::isfinite(scaled_double(v))) return HUGE_VAL;
        try { Q.push_back(to_exact(v)); } catch (const unreadable_value&) { return HUGE_VAL; }
    }
    double sum = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i; j < n; ++j) {
            exact m(i == j ? -1.0 : 0.0);
            for (int r = 0; r < n; ++r) m += Q[r * n + i] * Q[r * n + j];
            const double d = double(m);
            sum += (i == j ? 1.0 : 2.0) * d * d;
        }
    return std::sqrt(sum);
}

/// | ||y|| - ||y0|| | / ||y0|| for y = (u, w) and y0 = (u0, w0), exactly: how
/// far repeated application has stretched or shrunk y0. Measured against y0's
/// own length, so a start vector rounded off the unit circle adds nothing.
template <typename T>
double relative_length_change(const T& u, const T& w, const T& u0, const T& w0) {
    for (const T* v : {&u, &w})
        if (!std::isfinite(scaled_double(*v))) return HUGE_VAL;
    exact U, W, U0, W0;
    try {
        U = to_exact(u); W = to_exact(w); U0 = to_exact(u0); W0 = to_exact(w0);
    } catch (const unreadable_value&) {
        return HUGE_VAL;
    }
    exact d = U * U;
    d += W * W;
    exact d0 = U0 * U0;
    d0 += W0 * W0;
    const double n = std::sqrt(double(d)), n0 = std::sqrt(double(d0));
    d -= d0;
    return abs_double(d) / (n0 * (n + n0));   // | ||y||^2 - ||y0||^2 | / (||y0|| (||y|| + ||y0||))
}

} // namespace sw::mp_numerics
