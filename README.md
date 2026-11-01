# mp-numerics

[![CMake](https://github.com/stillwater-sc/mp-numerics/actions/workflows/cmake.yml/badge.svg)](https://github.com/stillwater-sc/mp-numerics/actions/workflows/cmake.yml)

**Mixed-precision numerical methods.** mp-numerics composes two header-only
libraries — [MTL5](https://github.com/stillwater-sc/mtl5) for linear algebra and
[Universal](https://github.com/stillwater-sc/universal) for parameterized number
systems — to study classic numerical methods (quadrature, root finding,
interpolation, optimization, approximation) under custom arithmetic (half
precision, posits, ...).

MTL5 deliberately has **no dependency on Universal**: it is the general
linear-algebra layer. mp-numerics is the integration layer where numerical
methods meet Universal's number types.

## Build

```bash
# Dependencies (MTL5 + Universal) are pulled automatically via FetchContent.
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j

# Run the smoke test
ctest --test-dir build --output-on-failure

# Run the quadrature precision study
./build/applications/integration/simpson_precision/simpson_precision

# Run the benchmark (optional arg: subinterval count)
./build/benchmarks/integration/simpson_benchmark 1048576
```

Using local checkouts instead of fetching from GitHub:

```bash
cmake -B build \
  -DFETCHCONTENT_SOURCE_DIR_MTL5=../mtl5 \
  -DFETCHCONTENT_SOURCE_DIR_UNIVERSAL=../universal
```

## Layout

The repo is organized by numerical-method category: **integration**, **roots**,
**interpolation**, **optimization**, **approximation**, **linalg**.

```
include/sw/mp_numerics/
  integration.hpp                # composite Simpson / trapezoidal, generic over Real
  exact.hpp                      # exact read-out of binary number types (ereal expansions)
  rotations.hpp reflectors.hpp   # Givens rotations and Householder reflectors, four scalings
  operator_metrics.hpp           # exact error measures of a rotation or reflector
  givens_qr.hpp householder_qr.hpp  # Givens QR; Householder QR storing beta as dgeqr2 does
applications/
  integration/simpson_precision/ # quadrature error across number precisions
benchmarks/
  integration/                   # accuracy-per-unit-time across number systems
  linalg/givens_householder/     # one rotation or reflector at a time, per type (#13)
  linalg/decompositions/         # QR and least squares across number types (#14)
tests/
  integration/                   # composite Simpson smoke test across types
  linalg/                        # exact read-out; rotations, reflectors, QR
  roots/ interpolation/ optimization/ approximation/   # (no tests yet -- roadmap)
docs/roadmap.md                  # milestones and migration plan
```

## License

MIT — see [LICENSE](LICENSE).
