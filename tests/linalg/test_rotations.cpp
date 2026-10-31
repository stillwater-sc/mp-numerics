// test_rotations: the three rotation forms (rotations.hpp), MTL5's own
// rotation and reflector, and the exact measures of operator_metrics.hpp, on
// cases whose outcome is known. Returns non-zero on failure.
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <random>

#include <mtl/mat/dense2D.hpp>
#include <mtl/vec/dense_vector.hpp>
#include <mtl/operation/givens.hpp>
#include <mtl/operation/householder.hpp>

#include <sw/mp_numerics/operator_metrics.hpp>
#include <sw/mp_numerics/reflectors.hpp>
#include <sw/mp_numerics/rotations.hpp>

#include <universal/number/cfloat/cfloat.hpp>
#include <universal/number/posit/posit.hpp>
#include <universal/number/qd/qd.hpp>

namespace {

using namespace sw::mp_numerics;

int failures = 0;

void check(bool ok, const char* what) {
    if (!ok) { std::cerr << "FAIL: " << what << '\n'; ++failures; }
    else std::cout << "ok:   " << what << '\n';
}

template <typename T>
operator_error rotation_error(const plane_rotation<T>& g, const T& a, const T& b, int k) {
    T y0 = a, y1 = b;
    apply_rotation(g, y0, y1);
    return measure_rotation(a, b, k, g.c, g.s, g.r, y0, y1);
}

template <typename T>
plane_rotation<T> mtl5_rotation(const T& a, const T& b) {
    mtl::mat::dense2D<T> H(2, 1);
    H(0, 0) = a;
    H(1, 0) = b;
    mtl::vec::dense_vector<T> rhs(2), c(1), s(1);
    rhs(0) = T(0);
    rhs(1) = T(0);
    mtl::apply_givens_rotation(H, rhs, c, s, 0);
    return {c(0), s(0), H(0, 0)};
}

template <typename T>
operator_error reflector_error(const T& a, const T& b, int k) {
    mtl::vec::dense_vector<T> x(2);
    x(0) = a;
    x(1) = b;
    auto [v, tau] = mtl::householder(x);
    const T v1 = v(1), dot = a + v1 * b;
    return measure_reflector(a, b, k, tau, v1, T(a - tau * dot), T(b - tau * v1 * dot));
}

template <typename T>
operator_error lapack_reflector_error(const T& a, const T& b, int k) {
    mtl::vec::dense_vector<T> x(2);
    x(0) = a;
    x(1) = b;
    auto h = householder_lapack_pow2(x);
    const T v1 = h.v(1), dot = a + v1 * b;
    return measure_reflector(a, b, k, h.tau, v1, T(a - h.tau * dot), T(b - h.tau * v1 * dot));
}

template <typename T>
operator_error textbook_reflector_error(const T& a, const T& b, int k) {
    mtl::vec::dense_vector<T> x(2);
    x(0) = a;
    x(1) = b;
    auto h = mtl::householder(x, mtl::tag::textbook_reflector{});
    T y1 = b;
    if (h.vtv > T(0)) y1 = b - ((T(2) * (h.v(0) * a + h.v(1) * b)) / h.vtv) * h.v(1);
    return measure_textbook_reflector(a, b, k, h.v(0), h.v(1), h.vtv, h.alpha, y1, k);
}

bool same_error(const operator_error& x, const operator_error& y) {
    return x.finite == y.finite && x.identity == y.identity && x.orth == y.orth && x.r_err == y.r_err &&
           x.annih == y.annih && x.annih_fp == y.annih_fp;
}

/// Every operator gives the same measures in H as in IEEE type F, at every
/// scale from F's subnormals to its overflow: cfloat rounds as IEEE 754 does,
/// so its results in the study can be read as the hardware type's would be.
template <typename F, typename H>
bool matches_ieee(std::mt19937_64& g) {
    using L = std::numeric_limits<F>;
    bool same = true;
    for (int k = L::min_exponent - L::digits; k < L::max_exponent; k += 3)
        for (int i = 0; i < 4; ++i) {
            const double th = double(g() >> 11) * 0x1.0p-53 * 6.283185307179586;
            const F fa = scale_pow2(F(std::cos(th)), k), fb = scale_pow2(F(std::sin(th)), k);
            if (fa == F(0) && fb == F(0)) continue;
            const H ha(fa), hb(fb);
            same = same && same_error(rotation_error(givens_unscaled(fa, fb), fa, fb, k),
                                      rotation_error(givens_unscaled(ha, hb), ha, hb, k))
                        && same_error(rotation_error(givens_safe_scaled(fa, fb), fa, fb, k),
                                      rotation_error(givens_safe_scaled(ha, hb), ha, hb, k))
                        && same_error(rotation_error(givens_pow2_scaled(fa, fb), fa, fb, k),
                                      rotation_error(givens_pow2_scaled(ha, hb), ha, hb, k))
                        && same_error(reflector_error(fa, fb, k), reflector_error(ha, hb, k))
                        && same_error(lapack_reflector_error(fa, fb, k), lapack_reflector_error(ha, hb, k))
                        && same_error(textbook_reflector_error(fa, fb, k), textbook_reflector_error(ha, hb, k));
        }
    return same;
}

} // namespace

int main() {
    const double eps = std::numeric_limits<double>::epsilon();
    std::mt19937_64 g(13);
    auto uni = [&] { return double(g() >> 11) * 0x1.0p-53 * 2 - 1; };

    {   // MTL5's rotation is givens_unscaled, bit for bit.
        bool same = true;
        for (int i = 0; i < 200; ++i) {
            const double a = uni(), b = uni();
            const auto m = mtl5_rotation(a, b), u = givens_unscaled(a, b);
            same = same && m.c == u.c && m.s == u.s && m.r == u.r;
        }
        check(same, "apply_givens_rotation matches givens_unscaled");
    }
    {   // In double, every form is orthogonal and accurate to a few eps. A
        // reflector's ||H^T H - I|| = |g| |g - 2|, g = tau (1 + v1^2), carries
        // the errors of tau (~3u) and twice those of v1 (~4u each), so its
        // bound is ~22 eps; a rotation's c^2 + s^2 - 1 stays near 2 eps.
        double rot = 0, hh_orth = 0, hh_other = 0;
        for (int i = 0; i < 500; ++i) {
            const double a = uni(), b = uni();
            for (const auto& e : {rotation_error(givens_unscaled(a, b), a, b, 0),
                                  rotation_error(givens_max_scaled(a, b), a, b, 0),
                                  rotation_error(givens_safe_scaled(a, b), a, b, 0)})
                rot = std::max({rot, e.orth, e.r_err, e.annih, e.annih_fp});
            const auto h = reflector_error(a, b, 0);
            hh_orth = std::max(hh_orth, h.orth);
            hh_other = std::max({hh_other, h.r_err, h.annih, h.annih_fp});
        }
        std::cout << "      worst double errors: rotations " << rot / eps << " eps, reflector orthogonality "
                  << hh_orth / eps << " eps, reflector r and annihilation " << hh_other / eps << " eps\n";
        check(rot <= 4 * eps, "double rotations: within 4 eps");
        check(hh_orth <= 24 * eps, "double reflector: orthogonality within 24 eps");
        check(hh_other <= 4 * eps, "double reflector: r and annihilation within 4 eps");
    }
    {   // a*a overflows: the unscaled form returns r = inf, the scaled ones do not.
        const double a = 1e200, b = 3e199;
        const int k = 665;   // 2^665 ~ 1.3e200
        check(!rotation_error(givens_unscaled(a, b), a, b, k).finite, "double 1e200: unscaled r overflows");
        const auto m = rotation_error(givens_max_scaled(a, b), a, b, k);
        const auto s = rotation_error(givens_safe_scaled(a, b), a, b, k);
        check(m.finite && m.r_err <= 2 * eps && m.orth <= 4 * eps, "double 1e200: max-scaled is accurate");
        check(s.finite && s.r_err <= 2 * eps && s.orth <= 4 * eps, "double 1e200: safe-scaled is accurate");
        const auto h = reflector_error(a, b, k);
        check(h.finite && h.orth <= 8 * eps && h.r_err <= 4 * eps, "double 1e200: householder is accurate");
    }
    {   // a*a underflows to zero: the unscaled form returns the identity, with
        // r = 0 in place of ||x||.
        const double a = 3e-200, b = 1e-200;
        const int k = -662;
        const auto u = rotation_error(givens_unscaled(a, b), a, b, k);
        check(u.finite && std::abs(u.r_err - 1.0) < 1e-12, "double 3e-200: unscaled r underflows to 0");
        const auto s = rotation_error(givens_safe_scaled(a, b), a, b, k);
        check(s.finite && s.r_err <= 2 * eps, "double 3e-200: safe-scaled is accurate");
    }
    {   // A tail below half an ulp of a: MTL5's householder returns the identity.
        const double a = 1.0, b = 1e-17;
        const auto h = reflector_error(a, b, 0);
        check(h.identity && std::abs(h.annih - b) <= b * eps, "householder: a negligible tail gives the identity");
    }
    {   // posit<16,2> saturates instead of overflowing: at 2^30, a*a = 2^60 goes
        // to maxpos = 2^56 and the unscaled r is wrong, but finite. The measure
        // must call it gross. The max-scaled r is right to the 4 fraction bits
        // a posit<16,2> has at 2^30.
        using P = sw::universal::posit<16, 2>;
        const P a(std::ldexp(1.0, 30)), b(std::ldexp(1.0, 29));
        const auto u = rotation_error(givens_unscaled(a, b), a, b, 30);
        check(u.finite && u.r_err > 0.5, "posit<16,2> 2^30: unscaled r saturates (gross)");
        const auto m = rotation_error(givens_max_scaled(a, b), a, b, 30);
        check(m.finite && m.r_err <= std::ldexp(1.0, -5), "posit<16,2> 2^30: max-scaled r to the posit's 4 bits");
    }
    {   // cfloat in IEEE layout reproduces float and double through every operator.
        using sw::universal::cfloat;
        check((matches_ieee<float, cfloat<32, 8, std::uint32_t, true, false, false>>(g)),
              "cfloat<32,8> gives float's measures, subnormals to overflow");
        check((matches_ieee<double, cfloat<64, 11, std::uint32_t, true, false, false>>(g)),
              "cfloat<64,11> gives double's measures, subnormals to overflow");
    }
    {   // A value the measures cannot read counts as gross. In a saturating
        // cfloat<80,15> at ||x|| = 2^9246, the textbook reflector's v^T v
        // (2^18492) clamps to max() ~ 2^16384; read against 2^(2 * 9246) it
        // is 2^-2108, below double's range.
        using S80 = sw::universal::cfloat<80, 15, std::uint32_t, true, false, true>;
        const S80 a = scale_pow2(S80(0.8), 9246), b = scale_pow2(S80(0.6), 9246);
        const auto t = textbook_reflector_error(a, b, 9246);
        check(t.finite && t.orth == HUGE_VAL && t.annih == HUGE_VAL,
              "cfloat<80,15,sat> 2^9246: clamped v^T v makes the textbook reflector gross");
        // qd can hold an exact residual near 2^-1000: it reads out, and the
        // reflector for x = (1, 2^-205) is accurate.
        using sw::universal::qd;
        const qd qa(1.0), qb(std::ldexp(1.3, -205));
        const auto q = lapack_reflector_error(qa, qb, 0);
        const double eq = double(std::numeric_limits<qd>::epsilon());
        check(q.finite && q.orth <= 8 * eq && q.annih_fp <= 8 * eq, "qd (1, 1.3*2^-205): dlarfg form reads out, accurate");
    }
    {   // The findings of mp-numerics#13, pinned.
        using P32 = sw::universal::posit<32, 2>;
        using P64 = sw::universal::posit<64, 2>;
        const double e32 = double(std::numeric_limits<P32>::epsilon());
        const double e64 = double(std::numeric_limits<P64>::epsilon());

        // MTL5's reflector, with Parlett's sign, puts v1 ~ -2 x0 / x1 where a
        // posit has few bits; dlarfg's sign keeps |v1| <= 1. x0 = 1, so
        // MTL5's scaling (a division by 1) rounds nothing: the sign is the cause.
        const P32 a(1.0), b(std::ldexp(1.3, -14));
        const auto m = reflector_error(a, b, 0), l = lapack_reflector_error(a, b, 0);
        std::cout << "      posit<32,2> (1, 1.3*2^-14): MTL5 reflector " << m.orth / e32
                  << " eps, dlarfg form " << l.orth / e32 << " eps\n";
        check(!m.identity && m.orth > 50 * e32, "posit<32,2> small tail: MTL5's reflector loses accuracy");
        check(l.orth <= 8 * e32, "posit<32,2> small tail: dlarfg's form stays within 8 eps");
        const auto t = textbook_reflector_error(a, b, 0);
        check(t.orth <= 8 * e32, "posit<32,2> small tail: the textbook reflector stays within 8 eps");

        // dlartg's thresholds guard range, not precision: inside them, a*a of a
        // posit<64,2> at 2^100 lands where the posit has few bits. Scaling by an
        // exact power of two keeps the squares near 1.
        const P64 f(std::ldexp(0.8, 100)), g(std::ldexp(0.6, 100));
        const auto sa = rotation_error(givens_safe_scaled(f, g), f, g, 100);
        const auto pw = rotation_error(givens_pow2_scaled(f, g), f, g, 100);
        std::cout << "      posit<64,2> 2^100: safe-scaled " << sa.orth / e64 << " eps, pow2-scaled "
                  << pw.orth / e64 << " eps\n";
        check(sa.finite && sa.orth > 1e3 * e64, "posit<64,2> 2^100: safe-scaled loses accuracy");
        check(pw.finite && pw.orth <= 4 * e64, "posit<64,2> 2^100: pow2-scaled stays within 4 eps");

        // The textbook reflector sums unscaled squares, so it overflows where
        // the unscaled rotation does; dlarfg's form with 2^-e scaling does not.
        const float fa = std::ldexp(0.8f, 100), fb = std::ldexp(0.6f, 100);
        check(!textbook_reflector_error(fa, fb, 100).finite, "float 2^100: textbook reflector overflows");
        const auto lf = lapack_reflector_error(fa, fb, 100);
        check(lf.finite && lf.orth <= 8 * double(std::numeric_limits<float>::epsilon()),
              "float 2^100: dlarfg form with 2^-e scaling is accurate");
    }
    return failures == 0 ? 0 : 1;
}
