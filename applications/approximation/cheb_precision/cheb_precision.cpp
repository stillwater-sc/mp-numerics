// cheb_precision: Chebyshev-series approximation accuracy across number
// precisions. Approximates 1/(1+x^2) on [-1,1] and reports the max error vs the
// series degree, per number type. Migrated in spirit from Universal's
// mixedprecision/approximation (Chebyshev study).
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string>

#include <sw/mp_numerics/approximation.hpp>

#include <universal/number/posit/posit.hpp>
#include <universal/number/cfloat/cfloat.hpp>

namespace {
double fx(double x) { return 1.0 / (1.0 + x * x); }

template <typename Real>
double cheb_error(std::size_t N) {
    sw::mp_numerics::chebyshev_series<Real> P([](Real x) { return Real(1) / (Real(1) + x * x); }, Real(-1), Real(1), N);
    double err = 0.0;
    for (int i = 0; i <= 200; ++i) { double t = -1.0 + 2.0 * i / 200.0; err = std::max(err, std::abs(double(P(Real(t))) - fx(t))); }
    return err;
}

template <typename Real>
void report(const std::string& tag) {
    std::cout << "  " << std::left << std::setw(14) << tag << std::right;
    for (std::size_t N : { std::size_t(8), std::size_t(16), std::size_t(32), std::size_t(48) })
        std::cout << "  N=" << std::setw(2) << N << " err=" << std::setw(10) << std::scientific << std::setprecision(2) << cheb_error<Real>(N);
    std::cout << '\n';
}
} // namespace

int main() {
    using namespace sw::universal;
    std::cout << "Chebyshev-series approximation of 1/(1+x^2) on [-1,1] -- max error vs degree\n";
    report<double>("double");
    report<float>("float");
    report<posit<32, 2>>("posit<32,2>");
    report<posit<16, 2>>("posit<16,2>");
    report<cfloat<16, 5, std::uint16_t, true, false, false>>("cfloat<16,5>");
    return 0;
}
