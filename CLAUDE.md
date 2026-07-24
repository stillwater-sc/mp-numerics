# CLAUDE.md

Guidance for Claude Code (claude.ai/code) when working in this repository.

## Project Overview

mp-numerics is the **integration layer** for mixed-precision numerical methods.
It composes two header-only sister libraries:

- [MTL5](https://github.com/stillwater-sc/mtl5) — C++20 linear algebra.
- [Universal](https://github.com/stillwater-sc/universal) — parameterized number
  systems (`cfloat`, `posit`, ...).

**Architectural rule:** MTL5 is the general linear-algebra layer and MUST NOT
depend on Universal. All MTL5 + Universal coupling lives here in mp-numerics.

## Build Commands

```bash
# Dependencies are pulled automatically via FetchContent.
cmake -B build -DCMAKE_BUILD_TYPE=Release -Wno-dev
cmake --build build -j
ctest --test-dir build --output-on-failure

# Use local sister checkouts instead of fetching from GitHub:
cmake -B build -DFETCHCONTENT_SOURCE_DIR_MTL5=../mtl5 \
               -DFETCHCONTENT_SOURCE_DIR_UNIVERSAL=../universal
```

## Architecture

- Header-only composition under `include/sw/mp_numerics/`. Namespace:
  `sw::mp_numerics`.
- CMake: INTERFACE library `sw::mp_numerics` linking MTL5 + Universal. Options:
  `MPNUMERICS_BUILD_APPLICATIONS`, `MPNUMERICS_BUILD_TESTS`,
  `MPNUMERICS_BUILD_BENCHMARKS`.
- `applications/`, `benchmarks/`, and `tests/` are organized by
  numerical-method category: `integration/` (quadrature), `roots/` (root
  finding), `interpolation/`, `optimization/`, `approximation/`.
- `applications/` — demonstration programs (each its own CMakeLists).
- `benchmarks/` — accuracy/performance studies.
- `tests/` — lightweight self-checking executables (no external framework);
  register with `mpnumerics_add_test`.
- `docs/roadmap.md` — milestones and the Universal->mp-numerics migration plan.

## Conventions

- C++20, header-only. Match the sister repos (mtl5, mp-iterative, mp-ir) for
  style and CMake structure.
- Conventional Commits. Feature branches + PRs to `main`; CI must pass.
- Never commit build artifacts or downloaded data.
