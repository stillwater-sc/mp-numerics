// mp-numerics smoke test: Newton divided-difference interpolation reproduces a
// low-degree polynomial exactly (up to roundoff) in each number type. Returns
// non-zero on failure (no external framework).
#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>

#include <sw/mp_numerics/interpolation.hpp>

#include <universal/number/posit/posit.hpp>
#include <universal/number/cfloat/cfloat.hpp>

namespace {
template <typename Real>
bool ok(const char* tag, double tol) {
    // exact cubic p(x) = x^3 - 2x + 1 sampled at 4 nodes -> interpolant == p
    auto p = [](double x) { return x * x * x - 2.0 * x + 1.0; };
    std::vector<Real> xs = { Real(-1), Real(0), Real(1), Real(2) }, ys;
    for (auto& x : xs) ys.push_back(Real(p(double(x))));
    sw::mp_numerics::newton_interpolant<Real> P(xs, ys);
    double err = 0.0;
    for (double t : { -0.5, 0.25, 0.7, 1.5 }) err = std::max(err, std::abs(double(P(Real(t))) - p(t)));
    if (err > tol) { std::cerr << tag << " interp failed: err=" << err << " tol=" << tol << '\n'; return false; }
    std::cout << tag << " interp ok: err=" << err << '\n';
    return true;
}
} // namespace

int main() {
    using namespace sw::universal;
    int failures = 0;
    if (!ok<double>("double", 1e-12))                                      ++failures;
    if (!ok<float>("float", 1e-4))                                         ++failures;
    if (!ok<posit<32, 2>>("posit<32,2>", 1e-6))                            ++failures;
    if (!ok<cfloat<32, 8, std::uint32_t, true, false, false>>("cfloat<32,8>", 1e-4)) ++failures;
    if (failures == 0) std::cout << "mp-numerics interp smoke test passed\n";
    return failures == 0 ? 0 : 1;
}
