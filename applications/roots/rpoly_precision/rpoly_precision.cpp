// mp-numerics -- Jenkins-Traub polynomial root finding (rpoly) accuracy across
// number precisions.
//
// For each test polynomial (well-separated real roots, a clustered/ill-separated
// case, and one with a complex-conjugate pair), find ALL roots via
// mtl::/sw::mp_numerics::poly_roots in a range of working precisions and report
// the worst error against the known roots. Illustrates where the achievable root
// accuracy is bounded by the working precision vs by the algorithm itself
// (universal#1211).

#include <algorithm>
#include <complex>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <sw/mp_numerics/polyroots.hpp>

#include <universal/number/posit/posit.hpp>
#include <universal/number/cfloat/cfloat.hpp>

namespace {

// worst distance from each known root to the nearest computed root (complex plane)
template <typename Real>
double max_root_error(const std::vector<sw::universal::complex<Real>>& computed,
                      const std::vector<std::complex<double>>& known) {
    double worst = 0.0;
    for (const auto& k : known) {
        double best = 1e300;
        for (const auto& c : computed) {
            const double dr = double(c.real()) - k.real();
            const double di = double(c.imag()) - k.imag();
            best = std::min(best, std::sqrt(dr * dr + di * di));
        }
        worst = std::max(worst, best);
    }
    return worst;
}

template <typename Real>
void row(const char* tag, const mtl::vec::dense_vector<Real>& coeffs,
         const std::vector<std::complex<double>>& known) {
    auto res = sw::mp_numerics::poly_roots<Real>(coeffs);
    const double err = max_root_error<Real>(res.roots, known);
    std::cout << "    " << std::setw(22) << std::left << tag
              << (res.ok ? "ok  " : "FAIL ")
              << "max root err = " << std::scientific << std::setprecision(3) << err << '\n';
}

// run one polynomial across a sweep of working precisions
template <typename Builder>
void study(const std::string& name, Builder build,
           const std::vector<std::complex<double>>& known) {
    using namespace sw::universal;
    std::cout << name << '\n';
    row<float>                                      ("float",        build.template operator()<float>(),        known);
    row<double>                                     ("double",       build.template operator()<double>(),       known);
    row<posit<16, 2>>                               ("posit<16,2>",  build.template operator()<posit<16, 2>>(),  known);
    row<posit<32, 2>>                               ("posit<32,2>",  build.template operator()<posit<32, 2>>(),  known);
    row<posit<64, 2>>                               ("posit<64,2>",  build.template operator()<posit<64, 2>>(),  known);
    row<cfloat<32, 8, std::uint32_t, true, false, false>>("cfloat<32,8>", build.template operator()<cfloat<32, 8, std::uint32_t, true, false, false>>(), known);
    std::cout << '\n';
}

} // namespace

int main()
try {
    using namespace sw::mp_numerics;
    using std::complex;

    std::cout << "Jenkins-Traub (rpoly) root accuracy vs working precision\n"
              << "========================================================\n\n";

    // (x-6)(x+1.5)(x-8): well separated real roots
    study("(x-6)(x+1.5)(x-8)  -- well-separated real roots",
          []<typename Real>() { return poly_from_real_roots<Real>({ Real(6), Real(-1.5), Real(8) }); },
          { {6, 0}, {-1.5, 0}, {8, 0} });

    // (x-1)(x-1.01)(x-1.02): tightly clustered real roots (precision-sensitive)
    study("(x-1)(x-1.01)(x-1.02)  -- clustered real roots",
          []<typename Real>() { return poly_from_real_roots<Real>({ Real(1), Real(1.01), Real(1.02) }); },
          { {1.0, 0}, {1.01, 0}, {1.02, 0} });

    // (x-2)(x^2+1) = x^3 - 2x^2 + x - 2: one real root + a complex conjugate pair
    study("(x-2)(x^2+1)  -- real root + complex pair",
          []<typename Real>() { return mtl::vec::dense_vector<Real>{ Real(1), Real(-2), Real(1), Real(-2) }; },
          { {2, 0}, {0, 1}, {0, -1} });

    return EXIT_SUCCESS;
}
catch (const std::exception& e) {
    std::cerr << "exception: " << e.what() << '\n';
    return EXIT_FAILURE;
}
catch (...) {
    std::cerr << "unknown exception\n";
    return EXIT_FAILURE;
}
