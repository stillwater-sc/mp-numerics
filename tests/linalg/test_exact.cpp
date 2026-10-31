// test_exact: values read out of binary number types exactly (exact.hpp).
// Returns non-zero on failure (no external framework).
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>

#include <sw/mp_numerics/exact.hpp>

#include <universal/number/fixpnt/fixpnt.hpp>
#include <universal/number/posit/posit.hpp>

namespace {

int failures = 0;

void check(bool ok, const char* what) {
    if (!ok) { std::cerr << "FAIL: " << what << '\n'; ++failures; }
    else std::cout << "ok:   " << what << '\n';
}

/// e == d exactly, for a double d.
bool equals(sw::mp_numerics::exact e, double d) {
    e -= sw::mp_numerics::exact(d);
    return double(e) == 0.0;
}

} // namespace

int main() {
    using namespace sw::mp_numerics;
    using sw::universal::posit;

    {   // More bits than a double: long double 1 + 2^-62 (x87; exact anywhere).
        const long double x = 1.0L + std::ldexp(1.0L, -std::numeric_limits<long double>::digits + 2);
        exact e = to_exact(x);
        e -= exact(1.0);
        check(equals(e, double(std::ldexp(1.0L, -std::numeric_limits<long double>::digits + 2))),
              "long double 1 + 2^-(p-2) reads out exactly");
    }
    {   // posit<64,2> carries up to 59 fraction bits: 1 + 2^-58.
        posit<64, 2> p(1);
        p += posit<64, 2>(std::ldexp(1.0, -58));
        exact e = to_exact(p);
        e -= exact(1.0);
        check(equals(e, std::ldexp(1.0, -58)), "posit<64,2> 1 + 2^-58 reads out exactly");
        const posit<64, 2> third = posit<64, 2>(1) / posit<64, 2>(3);
        exact t = to_exact(third);
        t -= exact(double(third));
        check(double(t) != 0.0, "posit<64,2> 1/3 holds bits beyond double's");
    }
    {   // Scaling past double's range: long double max * 2^-16383 = 2 - 2^-(p-1).
        if constexpr (std::numeric_limits<long double>::max_exponent > 1024) {
            const long double m = std::numeric_limits<long double>::max();
            const int k = std::numeric_limits<long double>::max_exponent - 1;
            exact e = to_exact(m, k);
            e -= exact(2.0);
            check(equals(e, -std::ldexp(1.0, 1 - std::numeric_limits<long double>::digits)),
                  "long double max, scaled by 2^-16383, reads out exactly");
        }
        const float d = std::numeric_limits<float>::denorm_min();
        check(equals(to_exact(d, -149), 1.0), "float denorm_min scaled by 2^149 is 1");
    }
    {   // An exact residual: c = 1 - 2^-30, s = 2^-15 give c^2 + s^2 - 1 = -2^-30 + 2^-60.
        const double c = 1.0 - std::ldexp(1.0, -30), s = std::ldexp(1.0, -15);
        exact o = to_exact(c) * to_exact(c);
        o += to_exact(s) * to_exact(s);
        o -= exact(1.0);
        check(equals(o, -std::ldexp(1.0, -30) + std::ldexp(1.0, -60)), "c^2 + s^2 - 1 is formed exactly");
    }
    {   // scale_pow2 in a type whose range is narrower than the shift.
        using F = sw::universal::fixpnt<16, 8, sw::universal::Modulo, std::uint16_t>;
        check(double(scale_pow2(F(1), -8)) == std::ldexp(1.0, -8), "fixpnt<16,8>: 1 * 2^-8 by small steps");
        check(double(scale_pow2(F(std::ldexp(1.0, -8)), 14)) == 64.0, "fixpnt<16,8>: 2^-8 * 2^14 by small steps");
        check(equals(to_exact(F(96.5), 6), 96.5 / 64), "fixpnt<16,8>: 96.5 scaled by 2^-6");
    }
    return failures == 0 ? 0 : 1;
}
