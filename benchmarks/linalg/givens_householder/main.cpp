// givens_householder_study: mp-numerics#13. Measures one Givens rotation or one
// Householder reflector at a time, in each number type, exactly, and writes
//   base.csv scale.csv tail.csv repeat.csv   (each with a .sysinfo sidecar)
// to the directory given as the first argument (default: the current one).
//
//   givens_householder_study [outdir] [--quick] [--family f1,f2,...]
//
// The families are ieee, posit, fixpnt, cfloat, lns, takum, dd and qd; all
// of them by default.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>

#include "study.hpp"

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

namespace {

using namespace study;

std::string num(double x) {
    if (std::isnan(x)) return "nan";
    if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
    char buf[32];
    std::snprintf(buf, sizeof buf, "%.6g", x);
    return buf;
}

double quantile(std::vector<double> v, double q) {
    if (v.empty()) return std::nan("");
    const std::size_t i = std::min(v.size() - 1, std::size_t(q * double(v.size() - 1) + 0.5));
    std::nth_element(v.begin(), v.begin() + std::ptrdiff_t(i), v.end());
    return v[i];
}

double maximum(const std::vector<double>& v) {
    return v.empty() ? std::nan("") : *std::max_element(v.begin(), v.end());
}

std::string compiler() {
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

std::string arch() {
#if defined(__x86_64__) || defined(_M_X64)
    return "x86_64";
#elif defined(__aarch64__) || defined(_M_ARM64)
    return "aarch64";
#else
    return "unknown";
#endif
}

std::string os() {
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

void write_sysinfo(const std::filesystem::path& csv, const config& cfg, bool quick) {
    std::ofstream f(csv.string() + ".sysinfo");
    const std::time_t now = std::time(nullptr);
    char date[32];
    std::strftime(date, sizeof date, "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&now));
    f << "label=" << csv.stem().string() << '\n'
      << "study=mp-numerics#13 givens_householder\n"
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
      << "long_double_digits=" << std::numeric_limits<long double>::digits << '\n'
      << "quick=" << (quick ? 1 : 0) << '\n'
      << "base_samples=" << cfg.base_samples << '\n'
      << "scale_points=" << cfg.scale_points << '\n'
      << "scale_samples=" << cfg.scale_samples << '\n'
      << "tail_samples=" << cfg.tail_samples << '\n'
      << "power_operators=" << cfg.power_operators << '\n'
      << "power_vectors=" << cfg.power_vectors << '\n'
      << "power_kmax=" << cfg.power_kmax << '\n'
      << "product_trials=" << cfg.product_trials << '\n'
      << "product_n=" << cfg.product_n << '\n'
      << "product_kmax=" << cfg.product_kmax << '\n';
}

void write_rows(const std::filesystem::path& path, const std::vector<row>& rows) {
    std::ofstream f(path);
    f << "type,family,bits,variant,eps,operator,input,param,samples,nonfinite,identity,gross,"
         "orth_max,orth_p50,r_err_max,r_err_p50,annih_max,annih_p50,annih_fp_max,annih_fp_p50\n";
    for (const row& w : rows) {
        const summary& s = w.s;
        f << '"' << w.type.name << "\"," << w.type.family << ',' << w.type.bits << ',' << w.type.variant << ','
          << num(w.eps) << ',' << w.op << ',' << w.input << ',' << num(w.param) << ','
          << s.samples << ',' << s.nonfinite << ',' << s.identity << ',' << s.gross << ','
          << num(maximum(s.orth)) << ',' << num(quantile(s.orth, 0.5)) << ','
          << num(maximum(s.r_err)) << ',' << num(quantile(s.r_err, 0.5)) << ','
          << num(maximum(s.annih)) << ',' << num(quantile(s.annih, 0.5)) << ','
          << num(maximum(s.annih_fp)) << ',' << num(quantile(s.annih_fp, 0.5)) << '\n';
    }
}

void write_repeat(const std::filesystem::path& path, const std::vector<repeat_row>& rows) {
    std::ofstream f(path);
    f << "type,family,bits,variant,eps,operator,mode,k,trials,nonfinite,dev_max,dev_p50\n";
    for (const repeat_row& w : rows)
        f << '"' << w.type.name << "\"," << w.type.family << ',' << w.type.bits << ',' << w.type.variant << ','
          << num(w.eps) << ',' << w.op << ',' << w.mode << ',' << w.k << ',' << w.trials << ',' << w.nonfinite << ','
          << num(maximum(w.dev)) << ',' << num(quantile(w.dev, 0.5)) << '\n';
}

/// The base class at a glance: worst errors in units of each type's epsilon.
void print_summary(const results& r) {
    std::printf("%-22s %-24s %10s %10s %10s %8s\n", "type", "operator", "orth/eps", "r_err/eps", "annih/eps", "gross");
    for (const row& w : r.base)
        std::printf("%-22s %-24s %10.3g %10.3g %10.3g %8zu\n", w.type.name.c_str(), w.op.c_str(),
                    maximum(w.s.orth) / w.eps, maximum(w.s.r_err) / w.eps, maximum(w.s.annih) / w.eps, w.s.gross);
}

} // namespace

int main(int argc, char** argv) try {
    std::filesystem::path outdir = ".";
    bool quick = false;
    std::string families = "ieee,posit,fixpnt,cfloat,lns,takum,dd,qd";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--quick") quick = true;
        else if (a == "--family" && i + 1 < argc) families = argv[++i];
        else outdir = a;
    }
    config cfg;
    if (quick) {
        cfg.base_samples = 200; cfg.scale_points = 24; cfg.scale_samples = 4; cfg.tail_samples = 4;
        cfg.power_operators = 2; cfg.power_vectors = 4; cfg.power_kmax = 256;
        cfg.product_trials = 1; cfg.product_kmax = 64;
    }
    std::filesystem::create_directories(outdir);

    results r;
    const auto t0 = std::chrono::steady_clock::now();
    auto want = [&](const char* f) { return ("," + families + ",").find(std::string(",") + f + ",") != std::string::npos; };
    if (want("ieee"))   run_ieee(cfg, r);
    if (want("posit"))  run_posit(cfg, r);
    if (want("fixpnt")) run_fixpnt(cfg, r);
    if (want("cfloat")) { run_cfloat_narrow(cfg, r); run_cfloat_mid(cfg, r); run_cfloat_wide(cfg, r); }
    if (want("lns"))    run_lns(cfg, r);
    if (want("takum"))  run_takum(cfg, r);
    if (want("dd"))     run_dd(cfg, r);
    if (want("qd"))     run_qd(cfg, r);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    const std::pair<const char*, const std::vector<row>*> tables[] = {
        {"base.csv", &r.base}, {"scale.csv", &r.scale}, {"tail.csv", &r.tail}};
    for (const auto& [name, rows] : tables) {
        write_rows(outdir / name, *rows);
        write_sysinfo(outdir / name, cfg, quick);
    }
    write_repeat(outdir / "repeat.csv", r.repeat);
    write_sysinfo(outdir / "repeat.csv", cfg, quick);

    print_summary(r);
    std::printf("wrote %s/{base,scale,tail,repeat}.csv in %.1f s\n", outdir.string().c_str(), sec);
    return 0;
} catch (const std::exception& e) {
    std::cerr << "givens_householder_study: " << e.what() << '\n';
    return 1;
}
