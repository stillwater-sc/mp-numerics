// interp_precision: Newton interpolation accuracy across number precisions.
// Interpolates 1/(1+x^2) at Chebyshev nodes over [-1,1] and reports the max
// error on a fine grid, per number type and node count. Migrated in spirit from
// Universal's mixedprecision/interpolation.
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string>
#include <numbers>
#include <vector>

#include <sw/mp_numerics/interpolation.hpp>

#include <universal/number/posit/posit.hpp>
#include <universal/number/cfloat/cfloat.hpp>

namespace {
double fx(double x) { return 1.0 / (1.0 + x * x); }

template <typename Real>
double interp_error(std::size_t nodes) {
    std::vector<Real> xs, ys;
    for (std::size_t k = 0; k < nodes; ++k) {
        double xk = std::cos(std::numbers::pi * (2.0 * k + 1.0) / (2.0 * nodes));   // Chebyshev node in [-1,1]
        xs.push_back(Real(xk)); ys.push_back(Real(fx(xk)));
    }
    sw::mp_numerics::newton_interpolant<Real> P(xs, ys);
    double err = 0.0;
    for (int i = 0; i <= 200; ++i) { double t = -1.0 + 2.0 * i / 200.0; err = std::max(err, std::abs(double(P(Real(t))) - fx(t))); }
    return err;
}

template <typename Real>
void report(const std::string& tag) {
    std::cout << "  " << std::left << std::setw(14) << tag << std::right;
    for (std::size_t n : { std::size_t(5), std::size_t(9), std::size_t(13), std::size_t(17) })
        std::cout << "  n=" << std::setw(2) << n << " err=" << std::setw(10) << std::scientific << std::setprecision(2) << interp_error<Real>(n);
    std::cout << '\n';
}
} // namespace

int main() {
    using namespace sw::universal;
    std::cout << "Newton interpolation of 1/(1+x^2) at Chebyshev nodes -- max error vs node count\n";
    report<double>("double");
    report<float>("float");
    report<posit<32, 2>>("posit<32,2>");
    report<posit<16, 2>>("posit<16,2>");
    report<cfloat<16, 5, std::uint16_t, true, false, false>>("cfloat<16,5>");
    return 0;
}
