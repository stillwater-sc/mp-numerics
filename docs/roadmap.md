# mp-numerics roadmap

## Milestone 0: composition layer bootstrapped (done)

- CMake scaffold replicated from the mp-* template (mp-iterative / mp-ir):
  INTERFACE library `sw::mp_numerics`, find_package -> FetchContent fallback for
  MTL5 + Universal, config-package install, CI matrix (MSVC/GCC/Clang/AppleClang).
- Shared header `include/sw/mp_numerics/integration.hpp`: composite Simpson and
  trapezoidal quadrature, generic over the arithmetic type.
- Smoke test: composite Simpson integrates known functions (x^2, 1/(1+x^2)) in
  `double`, `float`, `posit<32,2>`, `cfloat<32,8>` to each type's accuracy floor.
- Demo application `simpson_precision` (quadrature error vs subinterval count
  across number types) and benchmark `simpson_benchmark` (accuracy per unit
  time).
- Repo organized by numerical-method category: `integration/`, `roots/`,
  `interpolation/`, `optimization/`, `approximation/`.
- `include/mtl/math/quire_accumulator.hpp`: the MTL5 `accumulator_traits` <-
  Universal quire bridge (the coupling that must not live in MTL5).

## Milestone 1: migrate the Universal numerics experiments

Migrate from Universal `mixedprecision/{roots,integration,interpolation,
optimization,approximation}`, re-expressed on Universal precisions (and MTL5
kernels where the method uses linear algebra, e.g. Chebyshev approximation):

- [x] `roots` -- scalar root finding (bisection/Newton/secant,
  `include/sw/mp_numerics/roots.hpp`) across precisions: `root_precision` app +
  `test_root_smoke`. The Jenkins-Traub polynomial root finder (`rpoly`) is a
  larger follow-up.
- `integration` -- extend the quadrature study (adaptive rules, Gauss).
- [x] `interpolation` -- Newton divided-difference interpolation
  (`include/sw/mp_numerics/interpolation.hpp`): `interp_precision` app +
  `test_interp_smoke`. High-degree interpolation amplifies low-precision
  roundoff (posit<16,2>/cfloat<16,5> degrade as node count grows).
- [x] `optimization` -- golden-section 1-D minimization
  (`include/sw/mp_numerics/optimization.hpp`): `optim_precision` app +
  `test_optim_smoke`. Minimizer accuracy floors at ~sqrt(eps) of the type.
- `approximation` -- Chebyshev / minimax approximation (uses MTL5 dense LA).

## Milestone 2: mixed-precision numerics studies

- Where can low precision live inside an iterative numerical method without
  stalling? (mirrors the mp-ir / mp-iterative residual-precision question.)
- Quire (exact accumulation) inside quadrature / polynomial evaluation via the
  `accumulator_traits` bridge, vs naive same-precision accumulation.
