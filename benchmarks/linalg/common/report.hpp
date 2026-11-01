#pragma once
// Shared by the linear-algebra studies: number formatting and summary
// statistics for their CSVs, and the provenance lines that open every
// .sysinfo sidecar. The commits and flags come from mpn_study_provenance()
// in benchmarks/linalg/CMakeLists.txt.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <ctime>
#include <limits>
#include <ostream>
#include <string>
#include <thread>
#include <vector>

#ifndef MPN_GIT_COMMIT
#define MPN_GIT_COMMIT "unknown"
#endif
#ifndef MTL5_GIT_COMMIT
#define MTL5_GIT_COMMIT "unknown"
#endif
#ifndef UNIVERSAL_GIT_COMMIT
#define UNIVERSAL_GIT_COMMIT "unknown"
#endif
#ifndef MPN_BUILD_TYPE
#define MPN_BUILD_TYPE "unknown"
#endif
#ifndef MPN_CXX_FLAGS
#define MPN_CXX_FLAGS ""
#endif

namespace report {

/// x with six significant digits; nan, inf and -inf spelled out.
inline std::string num(double x) {
    if (std::isnan(x)) return "nan";
    if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.6g", x);
    return buf;
}

/// The q-quantile of v by nearest rank; NaN for an empty v.
inline double quantile(std::vector<double> v, double q) {
    if (v.empty()) return std::nan("");
    const std::size_t i = std::min(v.size() - 1, std::size_t(q * double(v.size() - 1) + 0.5));
    std::nth_element(v.begin(), v.begin() + std::ptrdiff_t(i), v.end());
    return v[i];
}

inline double maximum(const std::vector<double>& v) {
    return v.empty() ? std::nan("") : *std::max_element(v.begin(), v.end());
}

inline std::string compiler() {
#if defined(__clang__)
    return std::string("Clang ") + __clang_version__;
#elif defined(__GNUC__)
    return std::string("GCC ") + __VERSION__;
#elif defined(_MSC_VER)
    return "MSVC " + std::to_string(_MSC_VER);
#else
    return "unknown";
#endif
}

inline std::string arch() {
#if defined(__x86_64__) || defined(_M_X64)
    return "x86_64";
#elif defined(__aarch64__) || defined(_M_ARM64)
    return "aarch64";
#else
    return "unknown";
#endif
}

inline std::string os() {
#if defined(__linux__)
    return "Linux";
#elif defined(__APPLE__)
    return "macOS";
#elif defined(_WIN32)
    return "Windows";
#else
    return "unknown";
#endif
}

/// The sidecar lines every study writes first: what produced the file, where,
/// and from which sources. The study appends its own sample counts.
inline void write_provenance(std::ostream& f, const std::string& label, const std::string& study) {
    const std::time_t now = std::time(nullptr);
    char date[32];
    std::strftime(date, sizeof date, "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&now));
    f << "label=" << label << '\n'
      << "study=" << study << '\n'
      << "generated_utc=" << date << '\n'
      << "cpu_arch=" << arch() << '\n'
      << "cpu_logical_cores=" << std::thread::hardware_concurrency() << '\n'
      << "os_name=" << os() << '\n'
      << "compiler=" << compiler() << '\n'
      << "cpp_standard=" << __cplusplus << '\n'
      << "build_type=" << MPN_BUILD_TYPE << '\n'
      << "cxx_flags=" << MPN_CXX_FLAGS << '\n'
      << "git_commit=" << MPN_GIT_COMMIT << '\n'
      << "mtl5_commit=" << MTL5_GIT_COMMIT << '\n'
      << "universal_commit=" << UNIVERSAL_GIT_COMMIT << '\n'
      << "long_double_digits=" << std::numeric_limits<long double>::digits << '\n';
}

} // namespace report
