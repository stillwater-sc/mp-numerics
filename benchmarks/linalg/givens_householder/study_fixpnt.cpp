// givens_householder study: fixed point, with the binary point in the middle,
// wrapping (modulo) and saturating on overflow.
#include "study.hpp"

#include <universal/number/fixpnt/fixpnt.hpp>

namespace study {

void run_fixpnt(const config& cfg, results& out) {
    using sw::universal::fixpnt;
    using sw::universal::Modulo;
    using sw::universal::Saturate;
    run_all<fixpnt<16, 8, Modulo, std::uint16_t>>({"fixpnt<16,8,mod>", "fixpnt", "modulo", 16}, cfg, out);
    run_all<fixpnt<16, 8, Saturate, std::uint16_t>>({"fixpnt<16,8,sat>", "fixpnt", "saturate", 16}, cfg, out);
    run_all<fixpnt<32, 16, Modulo, std::uint32_t>>({"fixpnt<32,16,mod>", "fixpnt", "modulo", 32}, cfg, out);
    run_all<fixpnt<32, 16, Saturate, std::uint32_t>>({"fixpnt<32,16,sat>", "fixpnt", "saturate", 32}, cfg, out);
}

} // namespace study
