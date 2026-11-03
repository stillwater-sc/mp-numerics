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

## Documentation site

The site at https://stillwater-sc.github.io/mp-numerics/ is an Astro/Starlight
project in `docs-site/`, built and deployed to GitHub Pages by
`.github/workflows/docs.yml` (a pull request only builds it). It follows the
MTL5 and Universal sites.

**`docs-site/src/content/docs/` is 100% GENERATED — never write or edit files
there.** `docs-site/sync-content.mjs` wipes and regenerates it, together with
`docs-site/public/img/`, on every `npm run build` / `npm run dev`. Content
comes from two places:

- `docs/` — the site pages. `docs/site/*.mdx` are copied verbatim (they carry
  their own frontmatter, which must include an `editUrl` pointing at the
  `.mdx`); every other `docs/**/*.md` is listed in `FILE_MAP`,
  and gets its title from its H1 and its links rewritten.
- A study's `README.md` under `benchmarks/` — the README is the page, so it
  stays the single source next to the study's code. List the study directory
  in `STUDY_MAP`. Its `figures/` are copied to `public/img/studies/<slug>/`,
  and any other relative link (results, sources, scripts) is sent to the file
  on GitHub. Write figure links as `figures/<name>.svg` so they render on
  GitHub too.

To add a page: write it in `docs/<section>/` with an H1 and no frontmatter,
add it to `FILE_MAP`, and add it to `sidebar` in `docs-site/astro.config.mjs`
if the sidebar does not autogenerate that section. To add a study: add its
directory to `STUDY_MAP` and its slug to the Studies sidebar.

```bash
cd docs-site && npm ci && npm run dev    # preview at http://localhost:4321/mp-numerics/
```

## Conventions

- C++20, header-only. Match the sister repos (mtl5, mp-iterative, mp-ir) for
  style and CMake structure.
- Conventional Commits. Feature branches + PRs to `main`; CI must pass.
- Never commit build artifacts or downloaded data.
