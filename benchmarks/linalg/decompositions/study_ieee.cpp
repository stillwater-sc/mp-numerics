// decompositions study: the hardware IEEE types.
#include "study.hpp"

namespace study {

void run_ieee(const problems& P, const config& cfg, results& out) {
    run_all<float>({"float", "ieee", "binary32", 32}, P, cfg, out);
    run_all<double>({"double", "ieee", "binary64", 64}, P, cfg, out);
    // x87 80-bit on x86-64, binary128 on aarch64 Linux, double on MSVC and
    // Apple: the bits column records which one ran.
    constexpr int ld_bits = std::numeric_limits<long double>::digits == 64 ? 80
                          : std::numeric_limits<long double>::digits == 113 ? 128 : 64;
    run_all<long double>({"long double", "ieee", "long double", ld_bits}, P, cfg, out);
}

} // namespace study
