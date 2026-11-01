#pragma once
// givens_householder study (mp-numerics#13): one rotation or reflector at a
// time, in each number type, measured exactly. See README.md for the method.
//
// Seven operators, all built for x = (a, b):
//   givens_mtl5           mtl::apply_givens_rotation (operation/givens.hpp)
//   givens_max_scaled     sw::mp_numerics::givens_max_scaled (lapy2's scaling)
//   givens_safe_scaled    sw::mp_numerics::givens_safe_scaled (Anderson / dlartg)
//   givens_pow2_scaled    sw::mp_numerics::givens_pow2_scaled (exact 2^-e scaling)
//   householder_mtl5      mtl::householder (operation/householder.hpp), v(0) = 1
//   householder_textbook  mtl::householder(x, tag::textbook_reflector), v = x - alpha e1
//   householder_lapack_pow2 sw::mp_numerics::householder_lapack_pow2: dlarfg's sign
//                         and v(0) = 1, on x scaled by an exact 2^-e (the remedy tested)
//
// Four input classes:
//   base    random directions, ||x|| = 1
//   scale   random directions, ||x|| = 2^k over the type's whole range
//   tail    x = (+-1, +-2^-j (1 + u)), j from 0 to past the type's precision
//   repeat  one operator applied k times to unit vectors ("power"), and the
//           product of k fresh operators on an 8 x 8 Q ("product")

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <numbers>
#include <random>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <mtl/mat/dense2D.hpp>
#include <mtl/vec/dense_vector.hpp>
#include <mtl/operation/givens.hpp>
#include <mtl/operation/householder.hpp>

#include <universal/number/dd/dd_fwd.hpp>

#include <sw/mp_numerics/exact.hpp>
#include <sw/mp_numerics/operator_metrics.hpp>
#include <sw/mp_numerics/reflectors.hpp>
#include <sw/mp_numerics/rotations.hpp>

#include "common/random.hpp"

namespace study {

using sw::mp_numerics::operator_error;

struct type_desc {
    std::string name;     // e.g. posit<16,2>
    std::string family;   // ieee, posit, fixpnt, ...
    std::string variant;  // e.g. es2, modulo, saturate
    int bits;
    std::string seed = {};   // seeds the inputs; empty: the name. Variants of
                             // one format share it, so they see the same inputs.

    std::string seed_name() const { return seed.empty() ? name : seed; }
};

enum class op_kind {
    givens_mtl5, givens_max_scaled, givens_safe_scaled, givens_pow2_scaled,
    householder_mtl5, householder_textbook, householder_lapack_pow2
};
inline constexpr op_kind all_ops[] = {
    op_kind::givens_mtl5, op_kind::givens_max_scaled, op_kind::givens_safe_scaled,
    op_kind::givens_pow2_scaled, op_kind::householder_mtl5, op_kind::householder_textbook,
    op_kind::householder_lapack_pow2};

inline const char* op_name(op_kind op) {
    switch (op) {
    case op_kind::givens_mtl5:        return "givens_mtl5";
    case op_kind::givens_max_scaled:  return "givens_max_scaled";
    case op_kind::givens_safe_scaled: return "givens_safe_scaled";
    case op_kind::givens_pow2_scaled: return "givens_pow2_scaled";
    case op_kind::householder_mtl5:   return "householder_mtl5";
    case op_kind::householder_textbook: return "householder_textbook";
    case op_kind::householder_lapack_pow2: return "householder_lapack_pow2";
    }
    return "?";
}

/// Sample counts. The defaults produce the committed results; tests use less.
struct config {
    int base_samples    = 10000;
    int scale_points    = 160;   // grid of exponents k across the range
    int scale_samples   = 32;    // directions per k
    int tail_samples    = 32;    // per j
    int power_operators = 8;
    int power_vectors   = 16;
    int power_kmax      = 4096;
    int product_trials  = 4;
    int product_n       = 8;
    int product_kmax    = 1024;
};

/// Errors of many samples of one (type, operator, input, parameter).
struct summary {
    std::size_t samples = 0, nonfinite = 0, identity = 0, gross = 0;
    std::vector<double> orth, r_err, annih, annih_fp;

    void add(const operator_error& e, double eps) {
        ++samples;
        if (!e.finite) { ++nonfinite; return; }
        if (e.identity) ++identity;
        // Gross: the operator itself is wrong, beyond what rounding explains.
        // The 16 eps floor matters only for the 8-bit types, whose rounding
        // alone can put the errors above 1/2. r_err is left out, because at
        // the ends of a tapered or fixed-point range it mostly measures how
        // coarsely T represents ||x||, which no operator can improve.
        const double bound = std::max(0.5, 16 * eps);
        if (e.orth > bound || e.annih > bound) ++gross;
        orth.push_back(e.orth);
        r_err.push_back(e.r_err);
        annih.push_back(e.annih);
        annih_fp.push_back(e.annih_fp);
    }
};

struct row {
    type_desc type;
    double eps;
    std::string op, input;
    double param;
    summary s;
};

struct repeat_row {
    type_desc type;
    double eps;
    std::string op, mode;
    int k;
    std::size_t trials, nonfinite;
    std::vector<double> dev;   // one per trial, finite trials only
};

struct results {
    std::vector<row> base, scale, tail;
    std::vector<repeat_row> repeat;
};

// Reproducible randomness, shared with the other linear-algebra studies.
using report::seeded;

// ---- the operators, in T -----------------------------------------------------

/// A reflector with v(0) = 1, from MTL5's householder() or from dlarfg's form.
template <typename T>
std::pair<mtl::vec::dense_vector<T>, T> normalized(op_kind op, const mtl::vec::dense_vector<T>& x) {
    if (op == op_kind::householder_lapack_pow2) {
        auto h = sw::mp_numerics::householder_lapack_pow2(x);
        return {h.v, h.tau};
    }
    return mtl::householder(x);
}

/// The rotation forms that mp-numerics provides (not MTL5's own).
template <typename T>
sw::mp_numerics::plane_rotation<T> make_rotation(op_kind op, const T& a, const T& b) {
    using namespace sw::mp_numerics;
    if (op == op_kind::givens_max_scaled) return givens_max_scaled(a, b);
    if (op == op_kind::givens_pow2_scaled) return givens_pow2_scaled(a, b);
    return givens_safe_scaled(a, b);
}

template <typename T>
operator_error apply_and_measure(op_kind op, const T& a, const T& b, int k) {
    using namespace sw::mp_numerics;
    if (op == op_kind::householder_mtl5 || op == op_kind::householder_lapack_pow2) {
        mtl::vec::dense_vector<T> x(2);
        x(0) = a;
        x(1) = b;
        auto [v, tau] = normalized(op, x);
        const T v1 = v(1);
        // As qr_solve applies it: dot = v^T y, y(i) -= tau * v(i) * dot.
        const T dot = a + v1 * b;
        const T y0 = a - tau * dot;
        const T y1 = b - tau * v1 * dot;
        return measure_reflector(a, b, k, tau, v1, y0, y1);
    }
    if (op == op_kind::householder_textbook) {
        mtl::vec::dense_vector<T> x(2);
        x(0) = a;
        x(1) = b;
        auto h = mtl::householder(x, mtl::tag::textbook_reflector{});
        // As apply_householder_left applies it: f = 2 (v^T y) / vtv, y -= f v.
        T y1 = b;
        if (h.vtv > T(0)) {
            const T f = (T(2) * (h.v(0) * a + h.v(1) * b)) / h.vtv;
            y1 = b - f * h.v(1);
        }
        return measure_textbook_reflector(a, b, k, h.v(0), h.v(1), h.vtv, h.alpha, y1, k);
    }
    plane_rotation<T> g;
    if (op == op_kind::givens_mtl5) {
        mtl::mat::dense2D<T> H(2, 1);
        H(0, 0) = a;
        H(1, 0) = b;
        mtl::vec::dense_vector<T> rhs(2), c(1), s(1);
        rhs(0) = T(0);
        rhs(1) = T(0);
        mtl::apply_givens_rotation(H, rhs, c, s, 0);
        g = {c(0), s(0), H(0, 0)};
    } else {
        g = make_rotation(op, a, b);
    }
    T y0 = a, y1 = b;
    apply_rotation(g, y0, y1);
    return measure_rotation(a, b, k, g.c, g.s, g.r, y0, y1);
}

// ---- input classes -----------------------------------------------------------

template <typename T>
void run_base(const type_desc& t, const config& cfg, results& out) {
    const double eps = double(std::numeric_limits<T>::epsilon());
    for (op_kind op : all_ops) {
        seeded r(t.seed_name() + "/base");
        row w{t, eps, op_name(op), "base", 0, {}};
        for (int i = 0; i < cfg.base_samples; ++i) {
            const double th = r.angle();
            const T a(std::cos(th)), b(std::sin(th));
            if (a == T(0) && b == T(0)) continue;
            w.s.add(apply_and_measure(op, a, b, 0), eps);
        }
        out.base.push_back(std::move(w));
    }
}

/// The largest finite value of T. numeric_limits<dd>::max() is an infinity in
/// Universal, so dd takes double's, which it holds exactly.
template <typename T>
T largest() {
    if constexpr (std::is_same_v<T, sw::universal::dd>) return T(std::numeric_limits<double>::max());
    else return std::numeric_limits<T>::max();
}

/// The exponents k of the scale sweep: every k for a short range, otherwise
/// an even grid with the six exponents at each end added.
template <typename T>
std::vector<int> scale_grid(const config& cfg) {
    using L = std::numeric_limits<T>;
    using sw::mp_numerics::binary_exponent;
    const int lo = binary_exponent(L::denorm_min()), hi = binary_exponent(largest<T>());
    std::vector<int> ks;
    if (hi - lo <= cfg.scale_points) {
        for (int k = lo; k <= hi; ++k) ks.push_back(k);
        return ks;
    }
    for (int i = 0; i < cfg.scale_points; ++i) ks.push_back(lo + int(std::lround(double(hi - lo) * i / (cfg.scale_points - 1))));
    for (int d = 0; d < 6; ++d) { ks.push_back(lo + d); ks.push_back(hi - d); }
    std::sort(ks.begin(), ks.end());
    ks.erase(std::unique(ks.begin(), ks.end()), ks.end());
    return ks;
}

template <typename T>
void run_scale(const type_desc& t, const config& cfg, results& out) {
    using sw::mp_numerics::scale_pow2;
    const double eps = double(std::numeric_limits<T>::epsilon());
    const std::vector<int> ks = scale_grid<T>(cfg);
    for (op_kind op : all_ops) {
        seeded r(t.seed_name() + "/scale");
        for (int k : ks) {
            row w{t, eps, op_name(op), "scale", double(k), {}};
            for (int i = 0; i < cfg.scale_samples; ++i) {
                const double th = r.angle();
                const T a = scale_pow2(T(std::cos(th)), k), b = scale_pow2(T(std::sin(th)), k);
                if (a == T(0) && b == T(0)) continue;
                w.s.add(apply_and_measure(op, a, b, k), eps);
            }
            out.scale.push_back(std::move(w));
        }
    }
}

template <typename T>
void run_tail(const type_desc& t, const config& cfg, results& out) {
    const double eps = double(std::numeric_limits<T>::epsilon());
    const int jmax = -sw::mp_numerics::binary_exponent(std::numeric_limits<T>::epsilon()) + 8;
    for (op_kind op : all_ops) {
        seeded r(t.seed_name() + "/tail");
        for (int j = 0; j <= jmax; ++j) {
            row w{t, eps, op_name(op), "tail", double(j), {}};
            for (int i = 0; i < cfg.tail_samples; ++i) {
                const T a(r.uniform() < 0.5 ? -1.0 : 1.0);
                const double m = std::ldexp(1.0 + r.uniform(), -j);
                const T b(r.uniform() < 0.5 ? -m : m);
                w.s.add(apply_and_measure(op, a, b, 0), eps);
            }
            out.tail.push_back(std::move(w));
        }
    }
}

/// The operator built for (a, b), as a function applied in place to (u, w).
template <typename T>
std::function<void(T&, T&)> make_2x2(op_kind op, const T& a, const T& b) {
    using namespace sw::mp_numerics;
    if (op == op_kind::householder_mtl5 || op == op_kind::householder_lapack_pow2) {
        mtl::vec::dense_vector<T> x(2);
        x(0) = a;
        x(1) = b;
        auto [v, tau] = normalized(op, x);
        const T v1 = v(1), tv1 = tau * v1;
        return [tau, v1, tv1](T& u, T& w) {
            const T dot = u + v1 * w;
            u = u - tau * dot;
            w = w - tv1 * dot;
        };
    }
    if (op == op_kind::householder_textbook) {
        mtl::vec::dense_vector<T> x(2);
        x(0) = a;
        x(1) = b;
        auto h = mtl::householder(x, mtl::tag::textbook_reflector{});
        const T v0 = h.v(0), v1 = h.v(1), vtv = h.vtv;
        return [v0, v1, vtv](T& u, T& w) {
            if (!(vtv > T(0))) return;
            const T f = (T(2) * (v0 * u + v1 * w)) / vtv;
            u = u - f * v0;
            w = w - f * v1;
        };
    }
    plane_rotation<T> g;
    if (op == op_kind::givens_mtl5) {
        mtl::mat::dense2D<T> H(2, 1);
        H(0, 0) = a;
        H(1, 0) = b;
        mtl::vec::dense_vector<T> rhs(2), c(1), s(1);
        rhs(0) = T(0);
        rhs(1) = T(0);
        mtl::apply_givens_rotation(H, rhs, c, s, 0);
        g = {c(0), s(0), H(0, 0)};
    } else {
        g = make_rotation(op, a, b);
    }
    return [g](T& u, T& w) { apply_rotation(g, u, w); };
}

template <typename T>
void run_repeat(const type_desc& t, const config& cfg, results& out) {
    using namespace sw::mp_numerics;
    const double eps = double(std::numeric_limits<T>::epsilon());
    for (op_kind op : all_ops) {
        // power: one operator applied k times to unit vectors, measured
        // against their length as rounded into T.
        {
            seeded r(t.seed_name() + "/power");
            std::vector<std::vector<double>> dev;   // [checkpoint][operator]
            std::vector<std::size_t> nonfinite;
            std::vector<int> ks;
            for (int k = 1; k <= cfg.power_kmax; k *= 2) ks.push_back(k);
            dev.resize(ks.size());
            nonfinite.assign(ks.size(), 0);
            for (int o = 0; o < cfg.power_operators; ++o) {
                const double th = r.angle();
                auto apply = make_2x2<T>(op, T(std::cos(th)), T(std::sin(th)));
                std::vector<T> u(cfg.power_vectors), w(cfg.power_vectors);
                for (int i = 0; i < cfg.power_vectors; ++i) {
                    const double ph = r.angle();
                    u[i] = T(std::cos(ph));
                    w[i] = T(std::sin(ph));
                }
                const std::vector<T> u0 = u, w0 = w;
                std::size_t c = 0;
                for (int k = 1; k <= cfg.power_kmax; ++k) {
                    for (int i = 0; i < cfg.power_vectors; ++i) apply(u[i], w[i]);
                    if (k != ks[c]) continue;
                    double worst = 0;
                    for (int i = 0; i < cfg.power_vectors; ++i)
                        worst = std::max(worst, relative_length_change(u[i], w[i], u0[i], w0[i]));
                    if (std::isfinite(worst)) dev[c].push_back(worst);
                    else ++nonfinite[c];
                    ++c;
                }
            }
            for (std::size_t c = 0; c < ks.size(); ++c)
                out.repeat.push_back({t, eps, op_name(op), "power", ks[c],
                                      std::size_t(cfg.power_operators), nonfinite[c], dev[c]});
        }
        // product: Q <- G_k ... G_1 on an n x n Q, fresh operators each step.
        {
            seeded r(t.seed_name() + "/product");
            const int n = cfg.product_n;
            std::vector<int> ks;
            for (int k = 1; k <= cfg.product_kmax; k *= 2) ks.push_back(k);
            std::vector<std::vector<double>> dev(ks.size());
            std::vector<std::size_t> nonfinite(ks.size(), 0);
            for (int trial = 0; trial < cfg.product_trials; ++trial) {
                std::vector<T> Q(std::size_t(n) * n, T(0));
                for (int i = 0; i < n; ++i) Q[i * n + i] = T(1);
                std::size_t c = 0;
                for (int k = 1; k <= cfg.product_kmax; ++k) {
                    if (op == op_kind::householder_mtl5 || op == op_kind::householder_lapack_pow2) {
                        mtl::vec::dense_vector<T> x(n);
                        for (int i = 0; i < n; ++i) x(i) = T(2 * r.uniform() - 1);
                        auto [v, tau] = normalized(op, x);
                        if (tau != T(0)) {
                            for (int col = 0; col < n; ++col) {
                                T dot(0);
                                for (int i = 0; i < n; ++i) dot += v(i) * Q[i * n + col];
                                for (int i = 0; i < n; ++i) Q[i * n + col] -= tau * v(i) * dot;
                            }
                        }
                    } else if (op == op_kind::householder_textbook) {
                        mtl::vec::dense_vector<T> x(n);
                        for (int i = 0; i < n; ++i) x(i) = T(2 * r.uniform() - 1);
                        auto h = mtl::householder(x, mtl::tag::textbook_reflector{});
                        if (h.vtv > T(0)) {
                            for (int col = 0; col < n; ++col) {
                                T dot(0);
                                for (int i = 0; i < n; ++i) dot += h.v(i) * Q[i * n + col];
                                const T f = (T(2) * dot) / h.vtv;
                                for (int i = 0; i < n; ++i) Q[i * n + col] -= f * h.v(i);
                            }
                        }
                    } else {
                        const int i = r.below(n);
                        int j = r.below(n - 1);
                        if (j >= i) ++j;
                        const double th = r.angle();
                        auto apply = make_2x2<T>(op, T(std::cos(th)), T(std::sin(th)));
                        for (int col = 0; col < n; ++col) apply(Q[i * n + col], Q[j * n + col]);
                    }
                    if (k != ks[c]) continue;
                    const double d = orthogonality_frobenius(Q, n);
                    if (std::isfinite(d)) dev[c].push_back(d);
                    else ++nonfinite[c];
                    ++c;
                }
            }
            for (std::size_t c = 0; c < ks.size(); ++c)
                out.repeat.push_back({t, eps, op_name(op), "product", ks[c],
                                      std::size_t(cfg.product_trials), nonfinite[c], dev[c]});
        }
    }
}

template <typename T>
void run_all(const type_desc& t, const config& cfg, results& out) {
    run_base<T>(t, cfg, out);
    run_scale<T>(t, cfg, out);
    run_tail<T>(t, cfg, out);
    run_repeat<T>(t, cfg, out);
}

// One function per family, each in its own translation unit to keep the
// compile of any one file short.
void run_ieee(const config&, results&);
void run_posit(const config&, results&);
void run_fixpnt(const config&, results&);
void run_cfloat_narrow(const config&, results&);   // 8 to 32 bits
void run_cfloat_mid(const config&, results&);      // 40 to 64 bits
void run_cfloat_wide(const config&, results&);     // 80, 128, 256 bits
void run_lns(const config&, results&);
void run_takum(const config&, results&);
void run_dd(const config&, results&);
void run_qd(const config&, results&);

} // namespace study
