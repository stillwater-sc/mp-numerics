# Getting started

mp-numerics is header-only C++20. It composes two header-only sister
libraries:

- [MTL5](https://github.com/stillwater-sc/mtl5) for linear algebra;
- [Universal](https://github.com/stillwater-sc/universal) for parameterized
  number systems (`cfloat`, `posit`, `takum`, `lns`, `fixpnt`, `dd`, `qd`, ...).

MTL5 deliberately has no dependency on Universal. Everything that couples the
two lives here.

## Build

```bash
git clone https://github.com/stillwater-sc/mp-numerics.git
cd mp-numerics
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
ctest --test-dir build --output-on-failure
```

CMake looks for an installed Universal and an installed MTL5 of at least 5.22
through `find_package`. When it finds neither, it fetches both from `main` with
FetchContent.

| option | default | builds |
|---|---|---|
| `MPNUMERICS_BUILD_APPLICATIONS` | ON | the demonstration programs under `applications/` |
| `MPNUMERICS_BUILD_TESTS` | ON | the self-checking tests under `tests/` |
| `MPNUMERICS_BUILD_BENCHMARKS` | ON | the accuracy and performance studies under `benchmarks/` |

## Use local MTL5 and Universal checkouts

To measure a branch of MTL5 or Universal, or to reproduce a recorded run,
point CMake at the checkouts. The two `DISABLE` flags keep an installed copy
from taking precedence:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_DISABLE_FIND_PACKAGE_MTL5=ON -DCMAKE_DISABLE_FIND_PACKAGE_universal=ON \
  -DFETCHCONTENT_SOURCE_DIR_MTL5=../mtl5 \
  -DFETCHCONTENT_SOURCE_DIR_UNIVERSAL=../universal
```

The commits of all three repositories are read when CMake configures and are
written into each result's `.sysinfo` sidecar, so reconfigure before a run
whose results are kept.

## Run a study

Each study under `benchmarks/` is one executable, a `summarize.py` that prints
the tables in its README, and a `plot.py` that draws its figures. For the
Givens/Householder study:

```bash
cmake --build build --target givens_householder_study
./build/benchmarks/linalg/givens_householder/givens_householder_study results/
python3 benchmarks/linalg/givens_householder/summarize.py results/
python3 benchmarks/linalg/givens_householder/plot.py results/ figures/
```

`--quick` runs a small sample. `plot.py` needs matplotlib; `summarize.py`
needs only the standard library.

## Layout

```
include/sw/mp_numerics/   header-only composition (namespace sw::mp_numerics)
applications/             demonstration programs, by method category
benchmarks/               accuracy and performance studies
tests/                    self-checking executables
docs/                     this site's content
docs-site/                the Astro/Starlight site that publishes it
```
