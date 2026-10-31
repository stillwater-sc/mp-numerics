// givens_householder study: linear takums (base 2), with the specified three
// regime bits.
#include "study.hpp"

#include <universal/number/takum/takum.hpp>

namespace study {

void run_takum(const config& cfg, results& out) {
    using sw::universal::takum;
    run_all<takum<8>>({"takum<8>", "takum", "linear", 8}, cfg, out);
    run_all<takum<16>>({"takum<16>", "takum", "linear", 16}, cfg, out);
    run_all<takum<32>>({"takum<32>", "takum", "linear", 32}, cfg, out);
    run_all<takum<64>>({"takum<64>", "takum", "linear", 64}, cfg, out);
}

} // namespace study
