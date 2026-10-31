// givens_householder study: cfloat from 40 to 64 bits, all with double's
// 11-bit exponent.
#include "cfloat_variants.hpp"

namespace study {

void run_cfloat_mid(const config& cfg, results& out) {
    run_cfloat_variants<40, 11, std::uint32_t>(cfg, out);
    run_cfloat_variants<48, 11, std::uint32_t>(cfg, out);
    run_cfloat_variants<56, 11, std::uint32_t>(cfg, out);
    run_cfloat_variants<64, 11, std::uint32_t>(cfg, out);
}

} // namespace study
