// decompositions_study: mp-numerics#14, part 1. Factors the same qd-generated
// matrices in each number type by QR and solves least-squares problems by QR
// and by the normal equations, measuring the errors exactly, and writes
//   qr.csv lstsq.csv   (each with a .sysinfo sidecar)
// to the directory given as the first argument (default: the current one).
//
//   decompositions_study [outdir] [--quick] [--family f1,f2,...]
//
// The families are ieee, posit and fixpnt; all of them by default.
#include <chrono>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include "study.hpp"

#include "common/report.hpp"

namespace {

using namespace study;
using report::maximum;
using report::num;
using report::quantile;

void write_sysinfo(const std::filesystem::path& csv, const config& cfg, bool quick) {
    std::ofstream f(csv.string() + ".sysinfo");
    report::write_provenance(f, csv.stem().string(), "mp-numerics#14 decompositions part 1");
    auto list = [](const auto& v, auto&& show) {
        std::string s;
        for (const auto& e : v) s += (s.empty() ? "" : " ") + show(e);
        return s;
    };
    f << "quick=" << (quick ? 1 : 0) << '\n'
      << "samples=" << cfg.samples << '\n'
      << "kappa_decades=" << cfg.kappa_decades << '\n'
      << "kappa_step=" << cfg.kappa_step << '\n'
      << "kappa_past_eps=" << cfg.kappa_past_eps << '\n'
      << "qr_sizes=" << list(cfg.qr_sizes, [](std::size_t n) { return std::to_string(n); }) << '\n'
      << "ls_sizes=" << list(cfg.ls_sizes, [](const auto& p) {
             return std::to_string(p.first) + "x" + std::to_string(p.second);
         }) << '\n'
      << "kahan_thetas=" << list(cfg.kahan_thetas, [](double t) { return num(t); }) << '\n';
}

void write_qr(const std::filesystem::path& path, const std::vector<qr_row>& rows) {
    std::ofstream f(path);
    f << "type,family,bits,variant,eps,method,matrix,m,n,param,kappa,kappa_stored_p50,samples,nonfinite,thrown,"
         "q_identity,skipped,bwd_max,bwd_p50,bwd_col_max,bwd_col_p50,orth_max,orth_p50\n";
    for (const qr_row& w : rows)
        f << '"' << w.type.name << "\"," << w.type.family << ',' << w.type.bits << ',' << w.type.variant << ','
          << num(w.eps) << ',' << w.method << ',' << w.matrix << ',' << w.m << ',' << w.n << ',' << num(w.param)
          << ',' << num(w.kappa) << ',' << num(quantile(w.kappa_stored, 0.5)) << ',' << w.samples << ','
          << w.nonfinite << ',' << w.thrown << ',' << w.q_identity << ',' << w.skipped << ','
          << num(maximum(w.bwd)) << ',' << num(quantile(w.bwd, 0.5)) << ',' << num(maximum(w.bwd_col)) << ','
          << num(quantile(w.bwd_col, 0.5)) << ',' << num(maximum(w.orth)) << ',' << num(quantile(w.orth, 0.5))
          << '\n';
}

void write_ls(const std::filesystem::path& path, const std::vector<ls_row>& rows) {
    std::ofstream f(path);
    f << "type,family,bits,variant,eps,method,matrix,residual,m,n,param,kappa,kappa_stored_p50,samples,singular,"
         "nonfinite,thrown,breakdown,fwd_max,fwd_p50\n";
    for (const ls_row& w : rows)
        f << '"' << w.type.name << "\"," << w.type.family << ',' << w.type.bits << ',' << w.type.variant << ','
          << num(w.eps) << ',' << w.method << ',' << w.matrix << ',' << w.residual << ',' << w.m << ',' << w.n
          << ',' << num(w.param) << ',' << num(w.kappa) << ',' << num(quantile(w.kappa_stored, 0.5)) << ','
          << w.samples << ',' << w.singular << ',' << w.nonfinite << ',' << w.thrown << ',' << w.breakdown << ','
          << num(maximum(w.fwd)) << ',' << num(quantile(w.fwd, 0.5)) << '\n';
}

/// At a glance: QR of well-conditioned randsvd matrices at the largest size,
/// worst errors in units of each type's epsilon.
void print_summary(const results& r, const config& cfg) {
    const std::size_t n = cfg.qr_sizes.back();
    std::printf("randsvd, kappa = 1, n = %zu\n%-20s %-22s %10s %10s %10s\n", n, "type", "method", "bwd/eps",
                "bwdcol/eps", "orth/eps");
    for (const qr_row& w : r.qr)
        if (w.matrix == "randsvd" && w.param == 0 && w.n == n)
            std::printf("%-20s %-22s %10.3g %10.3g %10.3g\n", w.type.name.c_str(), w.method.c_str(),
                        maximum(w.bwd) / w.eps, maximum(w.bwd_col) / w.eps, maximum(w.orth) / w.eps);
}

} // namespace

int main(int argc, char** argv) try {
    std::filesystem::path outdir = ".";
    bool quick = false;
    std::string families = "ieee,posit,fixpnt";
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--quick") quick = true;
        else if (a == "--family" && i + 1 < argc) families = argv[++i];
        else outdir = a;
    }
    config cfg;
    if (quick) {
        cfg.samples = 2;
        cfg.kappa_step = 4;
        cfg.qr_sizes = {8};
        cfg.ls_sizes = {{16, 8}};
        cfg.kahan_thetas = {1.2, 0.4};
    }
    std::filesystem::create_directories(outdir);

    const auto t0 = std::chrono::steady_clock::now();
    const problems P = make_problems(cfg);
    const auto t1 = std::chrono::steady_clock::now();
    std::fprintf(stderr, "generated %zu QR and %zu least-squares groups in %.1f s\n", P.qr.size(), P.ls.size(),
                 std::chrono::duration<double>(t1 - t0).count());

    results r;
    auto want = [&](const char* f) { return ("," + families + ",").find(std::string(",") + f + ",") != std::string::npos; };
    if (want("ieee"))   run_ieee(P, cfg, r);
    if (want("posit"))  run_posit(P, cfg, r);
    if (want("fixpnt")) run_fixpnt(P, cfg, r);
    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    write_qr(outdir / "qr.csv", r.qr);
    write_sysinfo(outdir / "qr.csv", cfg, quick);
    write_ls(outdir / "lstsq.csv", r.ls);
    write_sysinfo(outdir / "lstsq.csv", cfg, quick);

    print_summary(r, cfg);
    std::printf("wrote %s/{qr,lstsq}.csv in %.1f s\n", outdir.string().c_str(), sec);
    return 0;
} catch (const std::exception& e) {
    std::cerr << "decompositions_study: " << e.what() << '\n';
    return 1;
}
