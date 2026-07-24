// simpson_precision: composite Simpson quadrature across number precisions.
//
// Integrates known functions with an increasing subinterval count and reports,
// per number type, the quadrature error against the exact value. Low-precision
// types stop improving once refinement error hits their rounding floor -- the
// mixed-precision quadrature story. Migrated in spirit from Universal's
// mixedprecision/integration (author's simpson study), re-expressed here on the
// mp-numerics composition layer.
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string>

#include <sw/mp_numerics/integration.hpp>

#include <universal/number/posit/posit.hpp>
#include <universal/number/cfloat/cfloat.hpp>

namespace {

// exact integral of 1/(1+x^2) over [0,1] is pi/4
const double kExact = std::atan(1.0);

template <typename Real>
void report(const std::string& tag) {
    using sw::mp_numerics::composite_simpson;
    auto f = [](Real x) { return Real(1) / (Real(1) + x * x); };
    std::cout << "  " << std::left << std::setw(14) << tag;
    for (std::size_t n : { std::size_t(8), std::size_t(32), std::size_t(128), std::size_t(512) }) {
        Real I = composite_simpson<Real>(f, Real(0), Real(1), n);
        double err = std::abs(double(I) - kExact);
        std::cout << "  n=" << std::setw(3) << std::right << n << " err=" << std::setw(10) << std::scientific << std::setprecision(2) << err;
    }
    std::cout << '\n';
}

} // namespace

int main() {
    using namespace sw::universal;
    std::cout << "composite Simpson quadrature of 1/(1+x^2) over [0,1] (exact = pi/4)\n";
    std::cout << "  number type      error vs subinterval count n\n";
    report<double>("double");
    report<float>("float");
    report<posit<32, 2>>("posit<32,2>");
    report<posit<16, 2>>("posit<16,2>");
    report<cfloat<32, 8, std::uint32_t, true, false, false>>("cfloat<32,8>");
    report<cfloat<16, 5, std::uint16_t, true, false, false>>("cfloat<16,5>");
    return 0;
}
