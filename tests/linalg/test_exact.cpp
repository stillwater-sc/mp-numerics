// test_exact: values read out of binary number types exactly (exact.hpp).
// Returns non-zero on failure (no external framework).
#include <cmath>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>

#include <sw/mp_numerics/exact.hpp>

#include <universal/number/cfloat/cfloat.hpp>
#include <universal/number/dd/dd.hpp>
#include <universal/number/fixpnt/fixpnt.hpp>
#include <universal/number/lns/lns.hpp>
#include <universal/number/posit/posit.hpp>
#include <universal/number/takum/takum.hpp>

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
    {   // binary_exponent beyond double's range, and where double(x) rounds up.
        using F128 = sw::universal::cfloat<128, 15, std::uint32_t, true, false, false>;
        using F256 = sw::universal::cfloat<256, 19, std::uint32_t, true, false, false>;
        check(binary_exponent(std::numeric_limits<F128>::max()) == 16383, "cfloat<128,15>: exponent of max is 16383");
        check(binary_exponent(std::numeric_limits<F128>::denorm_min()) == -16494,
              "cfloat<128,15>: exponent of denorm_min is -16494");
        check(binary_exponent(std::numeric_limits<F256>::denorm_min()) == -262378,
              "cfloat<256,19>: exponent of denorm_min is -262378");
        check(binary_exponent(F256(std::ldexp(1.5, 900)) * F256(std::ldexp(1.0, 900))) == 1800,
              "cfloat<256,19>: exponent of 1.5 * 2^1800 is 1800");
        posit<64, 2> p(1);
        p -= posit<64, 2>(std::ldexp(1.0, -58));
        check(double(p) == 1.0 && binary_exponent(p) == -1, "posit<64,2>: 1 - 2^-58 rounds to 1 in double, exponent -1");
        bool threw = false;
        try { (void)binary_exponent(std::numeric_limits<F128>::infinity()); } catch (const std::domain_error&) { threw = true; }
        check(threw, "binary_exponent rejects an infinity");
        using Tk = sw::universal::takum<8>;
        check(binary_exponent(std::numeric_limits<Tk>::min()) == -239, "takum<8>: exponent of min is -239");
    }
    {   // Wide cfloat values read out exactly, from beyond double's range too.
        using F256 = sw::universal::cfloat<256, 19, std::uint32_t, true, false, false>;
        F256 x(1);
        x += F256(std::ldexp(1.0, -60)) * F256(std::ldexp(1.0, -170));
        exact e = to_exact(x);
        e -= exact(1.0);
        check(equals(e, std::ldexp(1.0, -230)), "cfloat<256,19>: 1 + 2^-230 reads out exactly");
        exact m = to_exact(std::numeric_limits<F256>::max(), 262143);
        m -= exact(2.0);
        m += exact(std::ldexp(1.0, -236));
        check(double(m) == 0.0, "cfloat<256,19>: max scaled by 2^-262143 is 2 - 2^-236");
    }
    {   // dd scales by its own ldexp: its product turns DBL_MAX * 2^-512 into a NaN.
        using sw::universal::dd;
        const dd m(std::numeric_limits<double>::max());
        check(double(scale_pow2(m, -512)) == std::ldexp(std::numeric_limits<double>::max(), -512),
              "dd: DBL_MAX * 2^-512 is exact");
        check(binary_exponent(m) == 1023, "dd: exponent of DBL_MAX is 1023");
    }
    check(!is_binary_format<sw::universal::lns<16, 10, std::uint16_t>>::value, "lns is read through double");
    return failures == 0 ? 0 : 1;
}
