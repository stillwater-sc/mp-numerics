// optim_precision: 1-D golden-section minimization accuracy across precisions.
// Minimizes a unimodal function with a known minimizer and reports the
// minimizer accuracy per number type. Migrated in spirit from Universal's
// mixedprecision/optimization.
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <string>

#include <sw/mp_numerics/optimization.hpp>

#include <universal/number/posit/posit.hpp>
#include <universal/number/cfloat/cfloat.hpp>

namespace {
// f(x) = x^4 - 3x^3 + 2, unimodal on [0,3] with minimizer at x = 9/4 = 2.25
const double kMin = 2.25;

template <typename Real>
void report(const std::string& tag) {
    auto f = [](Real x) { return x * x * x * x - Real(3) * x * x * x + Real(2); };
    sw::mp_numerics::min_result r;
    Real xm = sw::mp_numerics::golden_section_min<Real>(f, Real(0), Real(3), Real(1e-10), 300, &r);
    std::cout << "  " << std::left << std::setw(14) << tag << std::right
              << "  minimizer err=" << std::setw(10) << std::scientific << std::setprecision(2) << std::abs(double(xm) - kMin)
              << "  iters=" << std::setw(3) << r.iterations
              << (r.converged ? "  converged" : "  (bracket floor)") << '\n';
}
} // namespace

int main() {
    using namespace sw::universal;
    std::cout << "golden-section minimization of x^4 - 3x^3 + 2 on [0,3] (min at x=2.25)\n";
    report<double>("double");
    report<float>("float");
    report<posit<32, 2>>("posit<32,2>");
    report<posit<16, 2>>("posit<16,2>");
    report<cfloat<16, 5, std::uint16_t, true, false, false>>("cfloat<16,5>");
    return 0;
}
