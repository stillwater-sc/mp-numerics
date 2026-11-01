#pragma once
// Reproducible randomness for the linear-algebra studies. std::mt19937_64 is
// specified bit for bit; std::hash and the distributions are not, so the seed
// is an FNV-1a hash of a name and uniforms come from the raw output.
#include <cstdint>
#include <numbers>
#include <random>
#include <string>

namespace report {

inline std::uint64_t fnv1a(const std::string& s) {
    std::uint64_t h = 1469598103934665603ull;
    for (unsigned char ch : s) { h ^= ch; h *= 1099511628211ull; }
    return h;
}

struct seeded {
    std::mt19937_64 g;
    explicit seeded(const std::string& seed) : g(fnv1a(seed)) {}
    double uniform() { return double(g() >> 11) * 0x1.0p-53; }   // [0, 1)
    double angle() { return 2 * std::numbers::pi * uniform(); }
    int below(int n) { return int(g() % std::uint64_t(n)); }
};

} // namespace report
