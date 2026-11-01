// decompositions study: posits, es = 2 as the 2022 standard fixes it, plus
// posit<16,1> for comparison with earlier results.
#include "study.hpp"

#include <universal/number/posit/posit.hpp>

namespace study {

void run_posit(const problems& P, const config& cfg, results& out) {
    using sw::universal::posit;
    run_all<posit<8, 2>>({"posit<8,2>", "posit", "es2", 8}, P, cfg, out);
    run_all<posit<16, 2>>({"posit<16,2>", "posit", "es2", 16}, P, cfg, out);
    run_all<posit<32, 2>>({"posit<32,2>", "posit", "es2", 32}, P, cfg, out);
    run_all<posit<64, 2>>({"posit<64,2>", "posit", "es2", 64}, P, cfg, out);
    run_all<posit<16, 1>>({"posit<16,1>", "posit", "es1", 16}, P, cfg, out);
}

} // namespace study
