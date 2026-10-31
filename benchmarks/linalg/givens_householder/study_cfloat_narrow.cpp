// givens_householder study: cfloat from 8 to 32 bits, including both FP8
// layouts (E4M3 and E5M2, without OCP's reuse of the top exponent).
#include "cfloat_variants.hpp"

namespace study {

void run_cfloat_narrow(const config& cfg, results& out) {
    run_cfloat_variants<8, 4, std::uint8_t>(cfg, out);
    run_cfloat_variants<8, 5, std::uint8_t>(cfg, out);
    run_cfloat_variants<16, 5, std::uint16_t>(cfg, out);
    run_cfloat_variants<24, 8, std::uint32_t>(cfg, out);
    run_cfloat_variants<32, 8, std::uint32_t>(cfg, out);
}

} // namespace study
