// mp-numerics smoke test: a Chebyshev series approximates a smooth function to
// the working precision's floor in each number type. Approximates 1/(1+x^2) on
// [-1,1] with degree 24. Returns non-zero on failure (no external framework).
#include <cmath>
#include <cstddef>
#include <iostream>

#include <sw/mp_numerics/approximation.hpp>

#include <universal/number/posit/posit.hpp>
#include <universal/number/cfloat/cfloat.hpp>

namespace {
double fx(double x) { return 1.0 / (1.0 + x * x); }

template <typename Real>
bool ok(const char* tag, double tol) {
    sw::mp_numerics::chebyshev_series<Real> P([](Real x) { return Real(1) / (Real(1) + x * x); }, Real(-1), Real(1), 24);
    double err = 0.0;
    for (int i = 0; i <= 100; ++i) { double t = -1.0 + 2.0 * i / 100.0; err = std::max(err, std::abs(double(P(Real(t))) - fx(t))); }
    if (err > tol) { std::cerr << tag << " cheb failed: err=" << err << " tol=" << tol << '\n'; return false; }
    std::cout << tag << " cheb ok: err=" << err << '\n';
    return true;
}
} // namespace

int main() {
    using namespace sw::universal;
    int failures = 0;
    if (!ok<double>("double", 1e-8))                                       ++failures;
    if (!ok<float>("float", 1e-3))                                         ++failures;
    if (!ok<posit<32, 2>>("posit<32,2>", 1e-4))                            ++failures;
    if (!ok<cfloat<32, 8, std::uint32_t, true, false, false>>("cfloat<32,8>", 1e-3)) ++failures;
    if (failures == 0) std::cout << "mp-numerics cheb smoke test passed\n";
    return failures == 0 ? 0 : 1;
}
