// root_precision: scalar root finding across number precisions.
//
// Finds a root with Newton and secant iterations and reports, per number type,
// the root accuracy against the exact value -- low-precision types stop short of
// the tolerance at their rounding floor. Migrated in spirit from Universal's
// mixedprecision/roots (secant study); the Jenkins-Traub polynomial root finder
// (rpoly) is a larger follow-up (see docs/roadmap.md).
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string>

#include <sw/mp_numerics/roots.hpp>

#include <universal/number/posit/posit.hpp>
#include <universal/number/cfloat/cfloat.hpp>

namespace {

// f(x) = x^3 - x - 2, real root ~ 1.5213797068045676
const double kRoot = 1.5213797068045676;

template <typename Real>
void report(const std::string& tag) {
    using namespace sw::mp_numerics;
    auto f  = [](Real x) { return x * x * x - x - Real(2); };
    auto df = [](Real x) { return Real(3) * x * x - Real(1); };

    root_result rn, rs;
    Real xn = newton<Real>(f, df, Real(2), Real(1e-12), 100, &rn);
    Real xs = secant<Real>(f, Real(1), Real(2), Real(1e-12), 100, &rs);

    std::cout << "  " << std::left << std::setw(14) << tag << std::right
              << "  newton: err=" << std::setw(10) << std::scientific << std::setprecision(2) << std::abs(double(xn) - kRoot)
              << " it=" << std::setw(3) << rn.iterations
              << "   secant: err=" << std::setw(10) << std::abs(double(xs) - kRoot)
              << " it=" << std::setw(3) << rs.iterations << '\n';
}

} // namespace

int main() {
    using namespace sw::universal;
    std::cout << "scalar root finding: x^3 - x - 2 = 0 (root ~ 1.52138), accuracy across precisions\n";
    report<double>("double");
    report<float>("float");
    report<posit<32, 2>>("posit<32,2>");
    report<posit<16, 2>>("posit<16,2>");
    report<cfloat<32, 8, std::uint32_t, true, false, false>>("cfloat<32,8>");
    report<cfloat<16, 5, std::uint16_t, true, false, false>>("cfloat<16,5>");
    return 0;
}
