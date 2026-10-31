// givens_householder study: logarithmic number systems, sized to the ranges of
// FP8 E4M3 (2^+-8), binary16 (2^+-16) and binary32 (2^+-128). Products and
// quotients are exact in lns; sums are rounded. lns values are not dyadic, so
// the measures read them through double (exact.hpp), which rounds 2^-29 below
// lns<32,23>'s own resolution.
#include "study.hpp"

#include <universal/number/lns/lns.hpp>

namespace study {

void run_lns(const config& cfg, results& out) {
    using sw::universal::lns;
    run_all<lns<8, 3, std::uint8_t>>({"lns<8,3>", "lns", "rbits3", 8}, cfg, out);
    run_all<lns<16, 10, std::uint16_t>>({"lns<16,10>", "lns", "rbits10", 16}, cfg, out);
    run_all<lns<32, 23, std::uint32_t>>({"lns<32,23>", "lns", "rbits23", 32}, cfg, out);
}

} // namespace study
