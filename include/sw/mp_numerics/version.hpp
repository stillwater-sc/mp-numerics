#pragma once
// mp-numerics -- mixed-precision numerical methods (MTL5 + Universal)
//
// Header-only composition layer. Shared numerical-method utilities live under
// sw::mp_numerics (roots, integration, interpolation, optimization,
// approximation). This header carries the version metadata.

namespace sw::mp_numerics {

inline constexpr int version_major = 0;
inline constexpr int version_minor = 1;
inline constexpr int version_patch = 0;

} // namespace sw::mp_numerics
