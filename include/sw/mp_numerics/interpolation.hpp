#pragma once
// mp-numerics -- polynomial interpolation (mixed precision).
//
// Newton divided-difference interpolation, generic over the arithmetic type so
// the same interpolant runs in double, float, posit, cfloat, ... and the study
// can observe how the working precision limits interpolation accuracy (and how
// the divided-difference table itself loses precision at high degree). Scalar
// arithmetic only.

#include <cstddef>
#include <vector>

namespace sw::mp_numerics {

/// Newton form of the interpolating polynomial through the points (x[i], y[i]).
/// Construction builds the divided-difference table; evaluation is Horner on the
/// Newton basis. O(n^2) build, O(n) evaluate.
template <typename Real>
class newton_interpolant {
public:
    newton_interpolant(const std::vector<Real>& xs, const std::vector<Real>& ys)
        : x_(xs), c_(ys) {
        const std::size_t n = x_.size();
        for (std::size_t j = 1; j < n; ++j)
            for (std::size_t i = n; i-- > j; )                 // i = n-1 .. j
                c_[i] = (c_[i] - c_[i - 1]) / (x_[i] - x_[i - j]);
    }

    Real operator()(Real t) const {
        const std::size_t n = c_.size();
        if (n == 0) return Real(0);
        Real p = c_[n - 1];
        for (std::size_t k = n - 1; k-- > 0; ) p = p * (t - x_[k]) + c_[k];
        return p;
    }

    std::size_t degree() const { return c_.empty() ? 0 : c_.size() - 1; }

private:
    std::vector<Real> x_;   // nodes
    std::vector<Real> c_;   // divided-difference coefficients
};

} // namespace sw::mp_numerics
