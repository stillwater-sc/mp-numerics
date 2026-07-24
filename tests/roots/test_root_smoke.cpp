// mp-numerics smoke test: bisection / Newton / secant find sqrt(2) (the root of
// x^2 - 2) in each number type, to the type's accuracy floor. Returns non-zero
// on failure (no external framework).
#include <cmath>
#include <cstddef>
#include <iostream>

#include <sw/mp_numerics/roots.hpp>

#include <universal/number/posit/posit.hpp>
#include <universal/number/cfloat/cfloat.hpp>

namespace {

template <typename Real>
bool ok(const char* tag, double tol) {
    using namespace sw::mp_numerics;
    auto f  = [](Real x) { return x * x - Real(2); };
    auto df = [](Real x) { return Real(2) * x; };
    const double sqrt2 = std::sqrt(2.0);

    Real rb = bisection<Real>(f, Real(1), Real(2), Real(1e-9), 200);
    Real rn = newton<Real>(f, df, Real(2), Real(1e-9), 100);
    Real rs = secant<Real>(f, Real(1), Real(2), Real(1e-9), 100);

    double eb = std::abs(double(rb) - sqrt2), en = std::abs(double(rn) - sqrt2), es = std::abs(double(rs) - sqrt2);
    bool good = (eb <= tol) && (en <= tol) && (es <= tol);
    if (!good) std::cerr << tag << " root failed: bisect=" << eb << " newton=" << en << " secant=" << es << " tol=" << tol << '\n';
    else       std::cout << tag << " roots ok: bisect=" << eb << " newton=" << en << " secant=" << es << '\n';
    return good;
}

} // namespace

int main() {
    using namespace sw::universal;
    int failures = 0;
    if (!ok<double>("double", 1e-8))                                       ++failures;
    if (!ok<float>("float", 1e-4))                                         ++failures;
    if (!ok<posit<32, 2>>("posit<32,2>", 1e-6))                            ++failures;
    if (!ok<cfloat<32, 8, std::uint32_t, true, false, false>>("cfloat<32,8>", 1e-4)) ++failures;
    if (failures == 0) std::cout << "mp-numerics root smoke test passed\n";
    return failures == 0 ? 0 : 1;
}
