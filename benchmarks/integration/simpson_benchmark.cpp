// simpson_benchmark: cost/accuracy of composite Simpson quadrature per number
// precision. Times a fixed-work integration and reports the achieved accuracy,
// so the accuracy-per-unit-time trade-off across number systems is visible.
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string>

#include <sw/mp_numerics/integration.hpp>

#include <universal/number/posit/posit.hpp>
#include <universal/number/cfloat/cfloat.hpp>

namespace {

const double kExact = std::atan(1.0);   // integral of 1/(1+x^2) over [0,1] = pi/4

template <typename Real>
void bench(const std::string& tag, std::size_t n) {
    using Clock = std::chrono::steady_clock;
    using sw::mp_numerics::composite_simpson;
    auto f = [](Real x) { return Real(1) / (Real(1) + x * x); };

    auto t0 = Clock::now();
    Real I = composite_simpson<Real>(f, Real(0), Real(1), n);
    double ms = std::chrono::duration<double, std::milli>(Clock::now() - t0).count();

    double err = std::abs(double(I) - kExact);
    std::cout << "  " << std::left << std::setw(14) << tag
              << "  n=" << n
              << "  time=" << std::setw(9) << std::fixed << std::setprecision(3) << ms << " ms"
              << "  err=" << std::scientific << std::setprecision(2) << err << '\n';
}

} // namespace

int main(int argc, char** argv) {
    using namespace sw::universal;
    std::size_t n = (argc > 1) ? static_cast<std::size_t>(std::stoul(argv[1])) : 1u << 20;
    std::cout << "composite Simpson benchmark -- 1/(1+x^2) over [0,1], n = " << n << " subintervals\n";
    std::cout << "  number type     work        time           accuracy\n";
    bench<double>("double", n);
    bench<float>("float", n);
    bench<posit<32, 2>>("posit<32,2>", n);
    bench<cfloat<32, 8, std::uint32_t, true, false, false>>("cfloat<32,8>", n);
    return 0;
}
