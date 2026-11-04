#pragma once
// givens_householder study: one cfloat size in its four variants. Exponent
// sizes follow IEEE 754 where it defines the size and truncate it elsewhere
// (cfloat<24,8>, cfloat<40..60,11>), so a size differs from its IEEE
// neighbours only in precision.
//   ieee    gradual underflow, overflow to infinity (IEEE 754 behaviour)
//   nosub   no subnormals: the smallest value is min(), as in many ML formats
//   sat     gradual underflow, overflow saturates to max() (as OCP's FP8 may)
//   maxexp  gradual underflow, overflow to infinity, and the all-ones exponent
//           holds finite values, as OCP's FP8 E4M3 reuses it: one more binade
//           at the top, with only its last encodings left for inf and NaN
#include <string>

#include "study.hpp"

#include <universal/number/cfloat/cfloat.hpp>

namespace study {

template <unsigned nbits, unsigned es, typename bt>
void run_cfloat_variants(const config& cfg, results& out) {
    using sw::universal::cfloat;
    const std::string n = "cfloat<" + std::to_string(nbits) + "," + std::to_string(es);
    // One seed for the four, so they see the same inputs.
    run_all<cfloat<nbits, es, bt, true, false, false>>({n + ">", "cfloat", "ieee", int(nbits), n + ">"}, cfg, out);
    run_all<cfloat<nbits, es, bt, false, false, false>>({n + ",nosub>", "cfloat", "nosub", int(nbits), n + ">"}, cfg, out);
    run_all<cfloat<nbits, es, bt, true, false, true>>({n + ",sat>", "cfloat", "sat", int(nbits), n + ">"}, cfg, out);
    run_all<cfloat<nbits, es, bt, true, true, false>>({n + ",maxexp>", "cfloat", "maxexp", int(nbits), n + ">"}, cfg, out);
}

} // namespace study
