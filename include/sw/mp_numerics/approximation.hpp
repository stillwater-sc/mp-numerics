#pragma once
// mp-numerics -- Chebyshev-series function approximation (mixed precision).
//
// Approximate f on [a, b] by a truncated Chebyshev series: sample f at the
// Chebyshev (Gauss) nodes, form the coefficients c_k, and evaluate by the
// Clenshaw recurrence. Generic over the arithmetic type, so the study can
// observe how the working precision limits the approximation -- the coefficient
// sums and the Clenshaw recurrence both accumulate in the working type. The node
// abscissae and cosine weights are computed in double (they are exact geometry,
// not part of the mixed-precision computation) and cast into the working type.

#include <cmath>
#include <cstddef>
#include <numbers>
#include <vector>

namespace sw::mp_numerics {

/// Truncated Chebyshev series of degree N approximating f on [a, b].
template <typename Real>
class chebyshev_series {
public:
    template <typename F>
    chebyshev_series(F f, Real a, Real b, std::size_t N) : a_(a), b_(b) {
        const std::size_t M = N + 1;                       // number of nodes
        std::vector<double> theta(M);
        std::vector<Real>   fx(M);
        for (std::size_t j = 0; j < M; ++j) {
            theta[j] = std::numbers::pi * (double(j) + 0.5) / double(M);
            const double xj = std::cos(theta[j]);                            // node in [-1,1]
            const Real node = Real(0.5) * (Real(xj) + Real(1)) * (b_ - a_) + a_;   // mapped to [a,b]
            fx[j] = f(node);
        }
        c_.assign(M, Real(0));
        for (std::size_t k = 0; k < M; ++k) {
            Real s(0);
            for (std::size_t j = 0; j < M; ++j) s = s + fx[j] * Real(std::cos(double(k) * theta[j]));
            c_[k] = (Real(2) / Real(static_cast<long long>(M))) * s;
        }
        c_[0] = c_[0] / Real(2);
    }

    /// Clenshaw evaluation at t in [a, b].
    Real operator()(Real t) const {
        const std::size_t M = c_.size();
        if (M == 0) return Real(0);
        const Real u = (Real(2) * t - a_ - b_) / (b_ - a_);   // map to [-1,1]
        Real bk1(0), bk2(0);
        for (std::size_t k = M - 1; k >= 1; --k) {
            const Real bk = Real(2) * u * bk1 - bk2 + c_[k];
            bk2 = bk1;
            bk1 = bk;
        }
        return u * bk1 - bk2 + c_[0];
    }

    std::size_t degree() const { return c_.empty() ? 0 : c_.size() - 1; }

private:
    std::vector<Real> c_;
    Real a_, b_;
};

} // namespace sw::mp_numerics
