// mp-numerics smoke test: the Jenkins-Traub polynomial root finder (rpoly)
// recovers the roots of known real polynomials, including a complex-conjugate
// pair. Tolerances are loose because this classic algorithm has a modest
// accuracy floor for well-separated roots that is independent of the working
// precision (see applications/roots/rpoly_precision). Returns non-zero on
// failure (no external framework).
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

#include <sw/mp_numerics/polyroots.hpp>
#include <universal/number/posit/posit.hpp>

namespace {

// worst distance from each known root to the nearest computed root
template <typename Real>
double worst_error(const std::vector<sw::universal::complex<Real>>& got,
                   const std::vector<std::complex<double>>& known) {
    double worst = 0.0;
    for (const auto& k : known) {
        double best = 1e300;
        for (const auto& c : got) {
            const double dr = double(c.real()) - k.real(), di = double(c.imag()) - k.imag();
            best = std::min(best, std::sqrt(dr * dr + di * di));
        }
        worst = std::max(worst, best);
    }
    return worst;
}

int failures = 0;
void check(bool cond, const std::string& what) {
    if (!cond) { ++failures; std::cout << "FAIL: " << what << '\n'; }
}

template <typename Real>
void run(const char* tag) {
    using namespace sw::mp_numerics;

    // (x-6)(x+1.5)(x-8): well-separated real roots
    {
        auto p = poly_from_real_roots<Real>({ Real(6), Real(-1.5), Real(8) });
        auto r = poly_roots<Real>(p);
        check(r.ok, std::string(tag) + ": separated cubic returned ok");
        check(r.roots.size() == 3, std::string(tag) + ": separated cubic has 3 roots");
        check(worst_error<Real>(r.roots, { {6,0},{-1.5,0},{8,0} }) < 0.1,
              std::string(tag) + ": separated cubic roots within 0.1");
    }
    // (x-2)(x^2+1): real root + complex conjugate pair
    {
        mtl::vec::dense_vector<Real> p{ Real(1), Real(-2), Real(1), Real(-2) };
        auto r = poly_roots<Real>(p);
        check(r.ok, std::string(tag) + ": complex-pair cubic returned ok");
        check(r.roots.size() == 3, std::string(tag) + ": complex-pair cubic has 3 roots");
        check(worst_error<Real>(r.roots, { {2,0},{0,1},{0,-1} }) < 0.1,
              std::string(tag) + ": complex-pair roots within 0.1");
    }
}

} // namespace

int main() {
    using namespace sw::universal;
    run<float>("float");
    run<double>("double");
    run<posit<32, 2>>("posit<32,2>");
    run<posit<64, 2>>("posit<64,2>");
    std::cout << (failures ? "rpoly smoke: FAILED\n" : "rpoly smoke: PASS\n");
    return failures ? 1 : 0;
}
