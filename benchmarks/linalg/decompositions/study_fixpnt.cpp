// decompositions study: fixed point, with the binary point in the middle,
// wrapping (modulo) and saturating on overflow.
//
// A fixpnt division by zero otherwise prints a message and returns 0 (modulo)
// or maxpos (saturate), which no later test can tell from a result. Thrown,
// it is counted in the CSVs' thrown column. Only this translation unit
// instantiates fixpnt, so the setting is the same wherever fixpnt is used.
#define FIXPNT_THROW_ARITHMETIC_EXCEPTION 1
#include "study.hpp"

#include <universal/number/fixpnt/fixpnt.hpp>

namespace study {

void run_fixpnt(const problems& P, const config& cfg, results& out) {
    using sw::universal::fixpnt;
    using sw::universal::Modulo;
    using sw::universal::Saturate;
    run_all<fixpnt<16, 8, Modulo, std::uint16_t>>({"fixpnt<16,8,mod>", "fixpnt", "modulo", 16}, P, cfg, out);
    run_all<fixpnt<16, 8, Saturate, std::uint16_t>>({"fixpnt<16,8,sat>", "fixpnt", "saturate", 16}, P, cfg, out);
    run_all<fixpnt<32, 16, Modulo, std::uint32_t>>({"fixpnt<32,16,mod>", "fixpnt", "modulo", 32}, P, cfg, out);
    run_all<fixpnt<32, 16, Saturate, std::uint32_t>>({"fixpnt<32,16,sat>", "fixpnt", "saturate", 32}, P, cfg, out);
}

} // namespace study
