// givens_householder study: cfloat wider than double, in range as well as
// precision: 80 and 128 bits with a 15-bit exponent, 256 with 19.
#include "cfloat_variants.hpp"

namespace study {

void run_cfloat_wide(const config& cfg, results& out) {
    run_cfloat_variants<80, 15, std::uint32_t>(cfg, out);
    run_cfloat_variants<128, 15, std::uint32_t>(cfg, out);
    run_cfloat_variants<256, 19, std::uint32_t>(cfg, out);
}

} // namespace study
