// mp-numerics smoke test: verify the MTL5 + Universal composition builds and
// that composite Simpson quadrature integrates known functions in each number
// type -- reaching the type's accuracy floor. Returns non-zero on failure (no
// external test framework, matching the mp-* repos' lightweight style).
#include <cmath>
#include <cstddef>
#include <iostream>

#include <sw/mp_numerics/integration.hpp>

#include <universal/number/posit/posit.hpp>
#include <universal/number/cfloat/cfloat.hpp>

namespace {

template <typename Real>
bool ok(const char* tag, double tol) {
    using sw::mp_numerics::composite_simpson;

    // (1) x^2 over [0,1] = 1/3. Simpson is exact for cubics, so only roundoff.
    Real i1 = composite_simpson<Real>([](Real x) { return x * x; }, Real(0), Real(1), 64);
    double e1 = std::abs(double(i1) - 1.0 / 3.0);

    // (2) 1/(1+x^2) over [0,1] = pi/4. O(h^4) truncation + roundoff.
    Real i2 = composite_simpson<Real>([](Real x) { return Real(1) / (Real(1) + x * x); }, Real(0), Real(1), 128);
    const double piover4 = std::atan(1.0);
    double e2 = std::abs(double(i2) - piover4);

    const bool good = (e1 <= tol) && (e2 <= tol);
    if (!good) std::cerr << tag << " Simpson failed: err(x^2)=" << e1 << " err(1/(1+x^2))=" << e2 << " tol=" << tol << '\n';
    else       std::cout << tag << " Simpson ok: err(x^2)=" << e1 << " err(1/(1+x^2))=" << e2 << '\n';
    return good;
}

} // namespace

int main() {
    using namespace sw::universal;
    int failures = 0;
    if (!ok<double>("double", 1e-12))                                   ++failures;
    if (!ok<float>("float", 1e-5))                                      ++failures;
    if (!ok<posit<32, 2>>("posit<32,2>", 1e-6))                         ++failures;
    if (!ok<cfloat<32, 8, std::uint32_t, true, false, false>>("cfloat<32,8>", 1e-5)) ++failures;
    if (failures == 0) std::cout << "mp-numerics smoke test passed\n";
    return failures == 0 ? 0 : 1;
}
