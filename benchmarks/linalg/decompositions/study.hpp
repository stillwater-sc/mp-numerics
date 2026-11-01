#pragma once
// decompositions study (mp-numerics#14, part 1): QR and least squares in each
// number type. See README.md for the method.
//
// Four QR methods, each returning an explicit Q (m x m) and R (m x n):
//   householder_mtl5      mtl::qr_factor, dlarfg's sign (MTL5 5.22.0)
//   householder_beta      sw::mp_numerics::qr_factor_stored_beta: the same
//                         reflectors, with beta stored on the diagonal as
//                         dgeqr2 does instead of recomputed (the remedy tested)
//   householder_textbook  mtl::qr_factor(A, tau, tag::textbook_reflector)
//   givens_mtl5           sw::mp_numerics::givens_qr, rotations from
//                         mtl::math::givens_rotation
//
// Four least-squares methods for min ||A x - b||:
//   qr_householder        mtl::qr_factor, then mtl::qr_solve
//   qr_householder_beta   qr_factor_stored_beta, then mtl::qr_solve
//   qr_givens             sw::mp_numerics::givens_qr_solve
//   normal_equations      A^T A and A^T b formed in T, then mtl::cholesky_factor
//                         and mtl::cholesky_solve
//
// Four matrix classes, generated in qd and rounded once into T, so every type
// factors the same matrices and everything is measured against the A the type
// actually stores:
//   randsvd   U S V^T, U and V random orthogonal, s_i = kappa^(-i/(n-1))
//   graded    U S: orthonormal columns scaled from 1 down to 1/kappa
//   kahan     U K, K = diag(s^i) (I - c N), s = sin(theta), c = cos(theta),
//             N strictly upper triangular ones; R is K up to signs
//   hilbert   1 / (i + j + 1)

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <mtl/mat/dense2D.hpp>
#include <mtl/operation/cholesky.hpp>
#include <mtl/operation/qr.hpp>
#include <mtl/operation/svd.hpp>
#include <mtl/vec/dense_vector.hpp>

#include <universal/number/qd/qd.hpp>

#include <sw/mp_numerics/exact.hpp>
#include <sw/mp_numerics/givens_qr.hpp>
#include <sw/mp_numerics/householder_qr.hpp>

#include "common/random.hpp"

namespace study {

using qd = sw::universal::qd;
template <typename T> using matrix = mtl::mat::dense2D<T>;
template <typename T> using vector = mtl::vec::dense_vector<T>;
using report::seeded;

struct type_desc {
    std::string name;     // e.g. posit<16,2>
    std::string family;   // ieee, posit, fixpnt
    std::string variant;  // e.g. es2, modulo, saturate
    int bits;
};

/// Sample counts and sizes. The defaults produce the committed results.
struct config {
    int samples = 8;          // random matrices per (class, size, parameter)
    int kappa_decades = 24;   // randsvd and graded: kappa = 10^0 .. 10^kappa_decades
    int kappa_step = 1;       // in decades
    int kappa_past_eps = 4;   // a type runs kappa up to 10^kappa_past_eps / eps
    std::vector<std::size_t> qr_sizes = {8, 16};
    std::vector<std::pair<std::size_t, std::size_t>> ls_sizes = {{16, 8}, {32, 16}};
    std::vector<double> kahan_thetas = {1.2, 0.8, 0.4, 0.2};
};

// ---- problems, in qd ---------------------------------------------------------

/// Matrices that share a class, a size and a parameter: one row of the CSVs.
struct qr_group {
    std::string cls;
    std::size_t m, n;
    double param;    // log10 kappa for randsvd and graded, theta for kahan, 0 for hilbert
    double kappa;    // of the qd matrices (the first one's, for kahan and hilbert)
    std::vector<matrix<qd>> A;
};

struct ls_group {
    std::string cls, residual;   // residual: "zero" (b = A x) or "large" (||r|| = ||A x||)
    std::size_t m, n;
    double param, kappa;
    std::vector<matrix<qd>> A;
    std::vector<vector<qd>> b;
};

struct problems {
    std::vector<qr_group> qr;
    std::vector<ls_group> ls;
};

inline qd uniform_pm1(seeded& r) { return qd(2 * r.uniform() - 1); }

/// A random orthogonal m x m matrix: Q of a matrix of uniforms, with R's
/// diagonal made positive so the distribution does not depend on the signs
/// qr_factor picks.
inline matrix<qd> random_orthogonal(std::size_t m, seeded& r) {
    matrix<qd> F(m, m);
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < m; ++j) F(i, j) = uniform_pm1(r);
    vector<qd> tau(m);
    mtl::qr_factor(F, tau);
    auto Q = mtl::qr_extract_Q(F, tau);
    auto R = mtl::qr_extract_R(F);
    mtl::qr_normalize_signs(Q, R);
    return Q;
}

/// s_i = 10^(-d i / (n - 1)), i = 0 .. n-1: from 1 down to 10^-d.
inline std::vector<qd> geometric(std::size_t n, double d) {
    std::vector<qd> s(n, qd(1));
    for (std::size_t i = 1; i < n; ++i)
        s[i] = sw::universal::pow(qd(10), qd(-d) * qd(double(i)) / qd(double(n - 1)));
    return s;
}

/// sigma_max / sigma_min of A (m >= n), in qd; infinity when A is singular.
inline double condition(const matrix<qd>& A) {
    auto [U, S, V] = mtl::svd(A);
    const std::size_t n = A.num_cols();
    qd lo = S(0, 0), hi = S(0, 0);
    for (std::size_t i = 1; i < n; ++i) {
        lo = std::min(lo, S(i, i));
        hi = std::max(hi, S(i, i));
    }
    return lo > qd(0) ? double(hi / lo) : HUGE_VAL;
}

inline matrix<qd> make_matrix(const std::string& cls, std::size_t m, std::size_t n, double param, seeded& r) {
    matrix<qd> A(m, n);
    if (cls == "hilbert") {
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < n; ++j) A(i, j) = qd(1) / qd(double(i + j + 1));
        return A;
    }
    const auto U = random_orthogonal(m, r);
    if (cls == "kahan") {
        const qd s = sw::universal::sin(qd(param)), c = sw::universal::cos(qd(param));
        matrix<qd> K(n, n);
        qd si(1);
        for (std::size_t i = 0; i < n; ++i) {
            for (std::size_t j = 0; j < n; ++j) K(i, j) = j < i ? qd(0) : (j == i ? si : -c * si);
            si *= s;
        }
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < n; ++j) {
                qd sum(0);
                for (std::size_t k = 0; k <= j; ++k) sum += U(i, k) * K(k, j);
                A(i, j) = sum;
            }
        return A;
    }
    const auto s = geometric(n, param);
    if (cls == "graded") {
        for (std::size_t i = 0; i < m; ++i)
            for (std::size_t j = 0; j < n; ++j) A(i, j) = U(i, j) * s[j];
        return A;
    }
    const auto V = random_orthogonal(n, r);   // randsvd
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = 0; j < n; ++j) {
            qd sum(0);
            for (std::size_t k = 0; k < n; ++k) sum += U(i, k) * s[k] * V(j, k);
            A(i, j) = sum;
        }
    return A;
}

/// b = A x for a random unit x, plus, for "large", an r orthogonal to A's
/// range with ||r|| = ||A x||, so the angle between b and range(A) is 45
/// degrees.
inline vector<qd> make_rhs(const matrix<qd>& A, const std::string& residual, seeded& r) {
    const std::size_t m = A.num_rows(), n = A.num_cols();
    vector<qd> x(n), b(m);
    qd nx(0);
    for (std::size_t j = 0; j < n; ++j) { x(j) = uniform_pm1(r); nx += x(j) * x(j); }
    nx = sw::universal::sqrt(nx);
    qd nb(0);
    for (std::size_t i = 0; i < m; ++i) {
        qd sum(0);
        for (std::size_t j = 0; j < n; ++j) sum += A(i, j) * x(j);
        b(i) = sum / nx;
        nb += b(i) * b(i);
    }
    if (residual == "zero") return b;
    matrix<qd> F(A);
    vector<qd> tau(n);
    mtl::qr_factor(F, tau);
    const auto Q = mtl::qr_extract_Q(F, tau);
    vector<qd> z(m), res(m);
    for (std::size_t i = 0; i < m; ++i) z(i) = uniform_pm1(r);
    qd nr(0);
    for (std::size_t i = 0; i < m; ++i) {   // res = (I - Q1 Q1^T) z
        res(i) = z(i);
        for (std::size_t k = 0; k < n; ++k) {
            qd p(0);
            for (std::size_t l = 0; l < m; ++l) p += Q(l, k) * z(l);
            res(i) -= Q(i, k) * p;
        }
        nr += res(i) * res(i);
    }
    const qd scale = sw::universal::sqrt(nb / nr);
    for (std::size_t i = 0; i < m; ++i) b(i) += scale * res(i);
    return b;
}

/// The parameters of the classes, for one size: (class, param, kappa when known).
inline std::vector<std::pair<std::string, double>> class_params(const config& cfg) {
    std::vector<std::pair<std::string, double>> p;
    for (const char* cls : {"randsvd", "graded"})
        for (int d = 0; d <= cfg.kappa_decades; d += cfg.kappa_step) p.emplace_back(cls, double(d));
    for (double th : cfg.kahan_thetas) p.emplace_back("kahan", th);
    p.emplace_back("hilbert", 0.0);
    return p;
}

inline problems make_problems(const config& cfg) {
    problems P;
    auto key = [](const std::string& what, const std::string& cls, std::size_t m, std::size_t n, double param,
                  int sample) {
        return what + "/" + cls + "/" + std::to_string(m) + "x" + std::to_string(n) + "/" + std::to_string(param) +
               "/" + std::to_string(sample);
    };
    for (std::size_t n : cfg.qr_sizes)
        for (const auto& [cls, param] : class_params(cfg)) {
            qr_group g{cls, n, n, param, 0, {}};
            const int samples = cls == "hilbert" ? 1 : cfg.samples;
            for (int i = 0; i < samples; ++i) {
                seeded r(key("qr", cls, n, n, param, i));
                g.A.push_back(make_matrix(cls, n, n, param, r));
            }
            g.kappa = (cls == "randsvd" || cls == "graded") ? std::pow(10.0, param) : condition(g.A.front());
            P.qr.push_back(std::move(g));
        }
    for (const auto& [m, n] : cfg.ls_sizes)
        for (const auto& [cls, param] : class_params(cfg))
            for (const char* residual : {"zero", "large"}) {
                ls_group g{cls, residual, m, n, param, 0, {}, {}};
                const int samples = cls == "hilbert" ? 1 : cfg.samples;
                for (int i = 0; i < samples; ++i) {
                    // The matrix does not depend on the residual, so both see the same A.
                    seeded r(key("ls", cls, m, n, param, i));
                    g.A.push_back(make_matrix(cls, m, n, param, r));
                    seeded rb(key(std::string("rhs-") + residual, cls, m, n, param, i));
                    g.b.push_back(make_rhs(g.A.back(), residual, rb));
                }
                g.kappa = (cls == "randsvd" || cls == "graded") ? std::pow(10.0, param) : condition(g.A.front());
                P.ls.push_back(std::move(g));
            }
    return P;
}

// ---- moving values between qd and T ------------------------------------------

/// x rounded into T: its leading two limbs, which hold 106 bits, summed in T.
/// That can round twice, so it is not always the nearest T, but nothing
/// depends on that: every error is measured against the value T stores.
template <typename T>
T from_qd(const qd& x) {
    if constexpr (std::is_same_v<T, qd>) return x;
    else {
        T t(x[0]);
        if (x[1] != 0.0) t = t + T(x[1]);
        return t;
    }
}

/// The exact value of a finite x of a binary type T, as a qd (exact.hpp's
/// peeling, into qd's 212 bits instead of an expansion).
template <typename T>
qd to_qd(const T& x) {
    T y = x;
    qd acc(0);
    for (int i = 0; i < 4 && y != T(0); ++i) {
        const double h = double(y);
        acc += qd(h);
        y = y - T(h);
    }
    return acc;
}

template <typename T>
bool finite(const T& x) {
    return std::isfinite(sw::mp_numerics::scaled_double(x));
}

template <typename T>
matrix<T> round_matrix(const matrix<qd>& A) {
    matrix<T> B(A.num_rows(), A.num_cols());
    for (std::size_t i = 0; i < A.num_rows(); ++i)
        for (std::size_t j = 0; j < A.num_cols(); ++j) B(i, j) = from_qd<T>(A(i, j));
    return B;
}

template <typename T>
matrix<qd> exact_matrix(const matrix<T>& A) {
    matrix<qd> B(A.num_rows(), A.num_cols());
    for (std::size_t i = 0; i < A.num_rows(); ++i)
        for (std::size_t j = 0; j < A.num_cols(); ++j) B(i, j) = to_qd(A(i, j));
    return B;
}

// ---- the methods, in T -------------------------------------------------------

enum class qr_method { householder_mtl5, householder_beta, householder_textbook, givens_mtl5 };
inline constexpr qr_method all_qr[] = {qr_method::householder_mtl5, qr_method::householder_beta,
                                       qr_method::householder_textbook, qr_method::givens_mtl5};

inline const char* name(qr_method q) {
    switch (q) {
    case qr_method::householder_mtl5:     return "householder_mtl5";
    case qr_method::householder_beta:     return "householder_beta";
    case qr_method::householder_textbook: return "householder_textbook";
    case qr_method::givens_mtl5:          return "givens_mtl5";
    }
    return "?";
}

enum class ls_method { qr_householder, qr_householder_beta, qr_givens, normal_equations };
inline constexpr ls_method all_ls[] = {ls_method::qr_householder, ls_method::qr_householder_beta,
                                       ls_method::qr_givens, ls_method::normal_equations};

inline const char* name(ls_method l) {
    switch (l) {
    case ls_method::qr_householder:   return "qr_householder";
    case ls_method::qr_householder_beta: return "qr_householder_beta";
    case ls_method::qr_givens:        return "qr_givens";
    case ls_method::normal_equations: return "normal_equations";
    }
    return "?";
}

template <typename T>
struct factors {
    matrix<T> Q, R;
    bool skipped = false;   // a reflector was left out (tau = 0) over a nonzero column
};

template <typename T>
factors<T> factor(qr_method method, const matrix<T>& A) {
    if (method == qr_method::givens_mtl5) {
        auto [Q, R] = sw::mp_numerics::givens_qr(A);
        return {std::move(Q), std::move(R), false};
    }
    const std::size_t m = A.num_rows(), n = A.num_cols();
    matrix<T> F(A);
    vector<T> tau(n);
    if (method == qr_method::householder_mtl5) mtl::qr_factor(F, tau);
    else if (method == qr_method::householder_beta) sw::mp_numerics::qr_factor_stored_beta(F, tau);
    else mtl::qr_factor(F, tau, mtl::tag::textbook_reflector{});
    bool skipped = false;
    for (std::size_t j = 0; j < std::min(m - 1, n); ++j) {
        if (tau(j) != T(0)) continue;
        for (std::size_t i = j + 1; i < m; ++i) skipped = skipped || F(i, j) != T(0);
    }
    return {mtl::qr_extract_Q(F, tau), mtl::qr_extract_R(F), skipped};
}

/// x, or empty when the method broke down (Cholesky found C not positive definite).
template <typename T>
std::pair<vector<T>, bool> solve(ls_method method, const matrix<T>& A, const vector<T>& b) {
    const std::size_t m = A.num_rows(), n = A.num_cols();
    vector<T> x(n);
    if (method == ls_method::qr_givens) return {sw::mp_numerics::givens_qr_solve(A, b), true};
    if (method == ls_method::qr_householder || method == ls_method::qr_householder_beta) {
        matrix<T> F(A);
        vector<T> tau(n);
        if (method == ls_method::qr_householder) mtl::qr_factor(F, tau);
        else sw::mp_numerics::qr_factor_stored_beta(F, tau);
        mtl::qr_solve(F, tau, x, b);
        return {x, true};
    }
    matrix<T> C(n, n);
    vector<T> d(n);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            T sum(0);
            for (std::size_t k = 0; k < m; ++k) sum += A(k, i) * A(k, j);
            C(i, j) = sum;
        }
        T sum(0);
        for (std::size_t k = 0; k < m; ++k) sum += A(k, i) * b(k);
        d(i) = sum;
    }
    if (mtl::cholesky_factor(C) != 0) return {x, false};
    mtl::cholesky_solve(C, x, d);
    return {x, true};
}

// ---- measures ----------------------------------------------------------------

/// The errors of one computed Q and R, each formed exactly and rounded once.
struct qr_error {
    double bwd;       ///< ||A - Q R||_F / ||A||_F
    double bwd_col;   ///< max_j ||(A - Q R) e_j|| / ||A e_j||: Householder's columnwise bound
    double orth;      ///< ||Q^T Q - I||_F
};

template <typename T>
qr_error measure(const matrix<T>& A, const factors<T>& f) {
    using sw::mp_numerics::exact;
    using sw::mp_numerics::to_exact;
    const std::size_t m = A.num_rows(), n = A.num_cols();
    std::vector<exact> XA(m * n), XQ(m * m), XR(m * n);
    try {
        for (std::size_t i = 0; i < m; ++i) {
            for (std::size_t j = 0; j < n; ++j) {
                XA[i * n + j] = to_exact(A(i, j));
                XR[i * n + j] = to_exact(f.R(i, j));
            }
            for (std::size_t j = 0; j < m; ++j) XQ[i * m + j] = to_exact(f.Q(i, j));
        }
    } catch (const sw::mp_numerics::unreadable_value&) {
        return {HUGE_VAL, HUGE_VAL, HUGE_VAL};
    }
    double total = 0, norm_a = 0, bwd_col = 0;
    for (std::size_t j = 0; j < n; ++j) {
        double col = 0, col_a = 0;
        for (std::size_t i = 0; i < m; ++i) {
            exact e = XA[i * n + j];
            for (std::size_t k = 0; k <= std::min(j, m - 1); ++k) e -= XQ[i * m + k] * XR[k * n + j];
            const double d = double(e), a = double(XA[i * n + j]);
            col += d * d;
            col_a += a * a;
        }
        total += col;
        norm_a += col_a;
        if (col_a > 0) bwd_col = std::max(bwd_col, std::sqrt(col / col_a));
        else if (col > 0) bwd_col = HUGE_VAL;
    }
    double orth = 0;
    for (std::size_t i = 0; i < m; ++i)
        for (std::size_t j = i; j < m; ++j) {
            exact e(i == j ? -1.0 : 0.0);
            for (std::size_t k = 0; k < m; ++k) e += XQ[k * m + i] * XQ[k * m + j];
            const double d = double(e);
            orth += (i == j ? 1.0 : 2.0) * d * d;
        }
    return {std::sqrt(total / norm_a), bwd_col, std::sqrt(orth)};
}

// ---- rows ----------------------------------------------------------------------

struct qr_row {
    type_desc type;
    double eps;
    std::string method, matrix;
    std::size_t m, n;
    double param, kappa;
    std::size_t samples = 0, nonfinite = 0, thrown = 0, q_identity = 0, skipped = 0;
    std::vector<double> kappa_stored = {}, bwd = {}, bwd_col = {}, orth = {};
};

struct ls_row {
    type_desc type;
    double eps;
    std::string method, matrix, residual;
    std::size_t m, n;
    double param, kappa;
    std::size_t samples = 0, singular = 0, nonfinite = 0, thrown = 0, breakdown = 0;
    std::vector<double> kappa_stored = {}, fwd = {};
};

struct results {
    std::vector<qr_row> qr;
    std::vector<ls_row> ls;
};

/// Whether a type runs a group: randsvd and graded stop kappa_past_eps
/// decades past 1/eps; kahan and hilbert always run.
inline bool in_reach(const std::string& cls, double param, double eps, const config& cfg) {
    if (cls != "randsvd" && cls != "graded") return true;
    return param <= std::log10(1 / eps) + cfg.kappa_past_eps;
}

template <typename T>
void run_qr(const type_desc& t, const problems& P, const config& cfg, results& out) {
    const double eps = double(std::numeric_limits<T>::epsilon());
    for (const qr_group& g : P.qr) {
        if (!in_reach(g.cls, g.param, eps, cfg)) continue;
        std::vector<qr_row> rows;
        for (qr_method q : all_qr) rows.push_back({t, eps, name(q), g.cls, g.m, g.n, g.param, g.kappa});
        for (const matrix<qd>& Aq : g.A) {
            const matrix<T> A = round_matrix<T>(Aq);
            const double ks = condition(exact_matrix(A));
            for (std::size_t q = 0; q < rows.size(); ++q) {
                qr_row& w = rows[q];
                ++w.samples;
                w.kappa_stored.push_back(ks);
                try {
                    const factors<T> f = factor(all_qr[q], A);
                    bool fin = true, ident = true;
                    for (std::size_t i = 0; i < g.m; ++i)
                        for (std::size_t j = 0; j < g.m; ++j) {
                            fin = fin && finite(f.Q(i, j));
                            ident = ident && f.Q(i, j) == (i == j ? T(1) : T(0));
                            if (j < g.n) fin = fin && finite(f.R(i, j));
                        }
                    if (f.skipped) ++w.skipped;
                    if (ident) ++w.q_identity;
                    if (!fin) { ++w.nonfinite; continue; }
                    const qr_error e = measure(A, f);
                    w.bwd.push_back(e.bwd);
                    w.bwd_col.push_back(e.bwd_col);
                    w.orth.push_back(e.orth);
                } catch (const std::exception&) {
                    ++w.thrown;
                }
            }
        }
        for (qr_row& w : rows) out.qr.push_back(std::move(w));
    }
}

template <typename T>
void run_ls(const type_desc& t, const problems& P, const config& cfg, results& out) {
    const double eps = double(std::numeric_limits<T>::epsilon());
    for (const ls_group& g : P.ls) {
        if (!in_reach(g.cls, g.param, eps, cfg)) continue;
        std::vector<ls_row> rows;
        for (ls_method l : all_ls)
            rows.push_back({t, eps, name(l), g.cls, g.residual, g.m, g.n, g.param, g.kappa});
        for (std::size_t s = 0; s < g.A.size(); ++s) {
            const matrix<T> A = round_matrix<T>(g.A[s]);
            vector<T> b(g.m);
            for (std::size_t i = 0; i < g.m; ++i) b(i) = from_qd<T>(g.b[s](i));
            // The reference: the least-squares solution of the problem T
            // stores, by Householder QR in qd.
            const matrix<qd> Aq = exact_matrix(A);
            const double ks = condition(Aq);
            vector<qd> bq(g.m), xref(g.n);
            for (std::size_t i = 0; i < g.m; ++i) bq(i) = to_qd(b(i));
            matrix<qd> F(Aq);
            vector<qd> tau(g.n);
            mtl::qr_factor(F, tau);
            qd rmin = abs(F(0, 0)), rmax = rmin;
            for (std::size_t k = 1; k < g.n; ++k) {
                rmin = std::min(rmin, abs(F(k, k)));
                rmax = std::max(rmax, abs(F(k, k)));
            }
            const bool singular = !(rmin > qd(1e-55) * rmax);
            if (!singular) mtl::qr_solve(F, tau, xref, bq);
            qd nref(0);
            for (std::size_t j = 0; j < g.n; ++j) nref += xref(j) * xref(j);
            for (std::size_t l = 0; l < rows.size(); ++l) {
                ls_row& w = rows[l];
                ++w.samples;
                w.kappa_stored.push_back(ks);
                if (singular) { ++w.singular; continue; }
                try {
                    const auto [x, ok] = solve(all_ls[l], A, b);
                    if (!ok) { ++w.breakdown; continue; }
                    bool fin = true;
                    for (std::size_t j = 0; j < g.n; ++j) fin = fin && finite(x(j));
                    if (!fin) { ++w.nonfinite; continue; }
                    qd d(0);
                    for (std::size_t j = 0; j < g.n; ++j) {
                        const qd e = to_qd(x(j)) - xref(j);
                        d += e * e;
                    }
                    w.fwd.push_back(double(sw::universal::sqrt(d / nref)));
                } catch (const std::exception&) {
                    ++w.thrown;
                }
            }
        }
        for (ls_row& w : rows) out.ls.push_back(std::move(w));
    }
}

template <typename T>
void run_all(const type_desc& t, const problems& P, const config& cfg, results& out) {
    run_qr<T>(t, P, cfg, out);
    run_ls<T>(t, P, cfg, out);
}

// One function per family, each in its own translation unit.
void run_ieee(const problems&, const config&, results&);
void run_posit(const problems&, const config&, results&);
void run_fixpnt(const problems&, const config&, results&);

} // namespace study
