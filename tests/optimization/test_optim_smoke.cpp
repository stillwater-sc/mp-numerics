// mp-numerics smoke test: golden-section search finds the minimizer of a
// unimodal function in each number type. Returns non-zero on failure.
//
// Note: golden section compares function values, so on a near-quadratic it can
// only locate the minimizer to ~sqrt(eps) (values flatten within sqrt(eps) of
// the minimum). The per-type tolerances reflect that floor, not full precision.
#include <cmath>
#include <cstddef>
#include <iostream>

#include <sw/mp_numerics/optimization.hpp>

#include <universal/number/posit/posit.hpp>
#include <universal/number/cfloat/cfloat.hpp>

namespace {
template <typename Real>
bool ok(const char* tag, double tol) {
    // f(x) = (x - 0.7)^2 + 1, minimizer at x = 0.7
    auto f = [](Real x) { return (x - Real(0.7)) * (x - Real(0.7)) + Real(1); };
    Real xm = sw::mp_numerics::golden_section_min<Real>(f, Real(0), Real(2), Real(1e-7), 200);
    double err = std::abs(double(xm) - 0.7);
    if (err > tol) { std::cerr << tag << " optim failed: err=" << err << " tol=" << tol << '\n'; return false; }
    std::cout << tag << " optim ok: minimizer err=" << err << '\n';
    return true;
}
} // namespace

int main() {
    using namespace sw::universal;
    int failures = 0;
    if (!ok<double>("double", 1e-6))                                       ++failures;
    if (!ok<float>("float", 2e-3))                                         ++failures;
    if (!ok<posit<32, 2>>("posit<32,2>", 2e-3))                            ++failures;
    if (!ok<cfloat<16, 5, std::uint16_t, true, false, false>>("cfloat<16,5>", 8e-2)) ++failures;
    if (failures == 0) std::cout << "mp-numerics optim smoke test passed\n";
    return failures == 0 ? 0 : 1;
}
