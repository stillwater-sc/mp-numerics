// givens_householder study: double-double and quad-double, whose range is
// double's and whose precision is two and four times its.
#include "study.hpp"

#include <universal/number/dd/dd.hpp>
#include <universal/number/qd/qd.hpp>

namespace study {

void run_dd(const config& cfg, results& out) {
    run_all<sw::universal::dd>({"dd", "dd", "double-double", 128}, cfg, out);
}

void run_qd(const config& cfg, results& out) {
    run_all<sw::universal::qd>({"qd", "qd", "quad-double", 256}, cfg, out);
}

} // namespace study
