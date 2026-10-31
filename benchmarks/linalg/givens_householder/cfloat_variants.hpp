#pragma once
// givens_householder study: one cfloat size in its three variants. Exponent
// sizes follow IEEE 754 where it defines the size and truncate it elsewhere
// (cfloat<24,8>, cfloat<40..56,11>), so a size differs from its IEEE
// neighbours only in precision.
//   ieee   gradual underflow, overflow to infinity (IEEE 754 behaviour)
//   nosub  no subnormals: the smallest value is min(), as in many ML formats
//   sat    gradual underflow, overflow saturates to max() (as OCP's FP8 may)
#include <string>

#include "study.hpp"

#include <universal/number/cfloat/cfloat.hpp>

namespace study {

template <unsigned nbits, unsigned es, typename bt>
void run_cfloat_variants(const config& cfg, results& out) {
    using sw::universal::cfloat;
    const std::string n = "cfloat<" + std::to_string(nbits) + "," + std::to_string(es);
    // One seed for the three, so they see the same inputs.
    run_all<cfloat<nbits, es, bt, true, false, false>>({n + ">", "cfloat", "ieee", int(nbits), n + ">"}, cfg, out);
    run_all<cfloat<nbits, es, bt, false, false, false>>({n + ",nosub>", "cfloat", "nosub", int(nbits), n + ">"}, cfg, out);
    run_all<cfloat<nbits, es, bt, true, false, true>>({n + ",sat>", "cfloat", "sat", int(nbits), n + ">"}, cfg, out);
}

} // namespace study
