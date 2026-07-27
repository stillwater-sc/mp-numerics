# Jenkins-Traub coefficient scaling and tapered-precision number types

`rpoly_ak1.hpp` is the canonical ak1 / Netlib TOMS-493 Jenkins-Traub real-polynomial
root finder, templatized over the arithmetic type `Real` so the same solver runs in
`double`, `float`, `posit`, `cfloat`, ... The single non-mechanical change from the
reference implementation is a **precision-aware coefficient-scaling gate**. This note
records why it is needed.

## The problem

The classic Jenkins-Traub algorithm begins by **unconditionally** rescaling the
polynomial coefficients by a power of two:

```
lo  = min() / epsilon();          // smallest "safe" magnitude for the type
sc  = lo / moduli_min;            // moduli_min = smallest |coefficient|
factor = 2^round(log2(sc));
for (i) p[i] *= factor;           // applied whenever a coefficient is "large enough"
```

The intent is numerical robustness: keep the coefficients away from the type's
overflow and underflow limits during the shift recurrences. For a **uniform-precision**
IEEE format this is harmless -- a power-of-two scale is exact and every representable
magnitude carries the same number of significant bits, so shifting the whole polynomial
up or down the exponent range costs nothing.

For a **tapered-precision** format it is catastrophic.

A posit does not carry uniform precision across its range: it has the most fraction bits
near magnitude 1 and progressively fewer toward its (very large) dynamic-range extremes.
Two consequences combine:

1. `min()` is astronomically small and `epsilon()` astronomically small for a posit,
   so `lo = min()/epsilon()` is a tiny number and `factor` comes out around `2^-92`
   for `posit<32,2>`.
2. Multiplying O(1) coefficients by `2^-92` moves them from the high-precision region
   near 1 deep into the low-precision tapered region.

The root finder then does all of its work on coefficients that have been stripped of
most of their significant bits. Accuracy collapses:

| coefficient scaling | `posit<32,2>` max root err, (x-6)(x+1.5)(x-8) |
|---------------------|-----------------------------------------------|
| unconditional (classic) | 2.5e-1 |
| gated (this header) | 1.2e-7 |

## The fix

Rescale **only when the coefficients actually approach the type's representable
extremes** -- which is the only situation the scaling was ever meant to protect
against:

```cpp
const Real hi_guard = sqrt(std::numeric_limits<Real>::max());
const Real lo_guard = sqrt(std::numeric_limits<Real>::min());
if ((moduli_max > hi_guard || moduli_min < lo_guard) && /* original trigger */) {
    // ... compute factor and apply ...
}
```

- **In-range coefficients** (the common case, including every polynomial in the
  precision study) are left untouched, so each type keeps its full precision near
  magnitude 1.
- **Extreme-magnitude polynomials** -- coefficients above `sqrt(max())` or below
  `sqrt(min())`, where unscaled evaluation would genuinely over/underflow -- still get
  the protective rescale, on any type.

The guard band is deliberately generous (`sqrt` of the extremes) so that scaling only
fires when it is actually needed to avoid overflow/underflow, never merely because a
coefficient happens to be "large" in absolute terms.

## Effect

With the gate, every number type reaches its precision floor instead of a fixed
scaling-induced floor:

| polynomial | double | posit<64,2> | posit<32,2> | float / cfloat<32,8> |
|------------|--------|-------------|-------------|----------------------|
| (x-6)(x+1.5)(x-8)     | 3.6e-15 | **0.0 (exact)** | 1.2e-7 | 1.9e-6 |
| (x-1)(x-1.01)(x-1.02) | 3.6e-12 | 9.5e-15 | 6.3e-4 | 6.3e-4 |
| (x-2)(x^2+1)          | 1.4e-17 | 2.3e-19 | 5.8e-10 | 7.5e-9 |

The precision study (`applications/roots/rpoly_precision`) now shows clean, monotonic
precision scaling -- higher working precision yields strictly better roots. This
relationship was previously masked: the earlier lower-fidelity port floored at ~1e-2
for well-separated `double` roots regardless of the type's precision.

## Takeaway

This is a general lesson for porting classic numerical kernels to tapered / non-uniform
number systems: **any preconditioning step that assumes uniform precision across the
exponent range** (power-of-two scaling, balancing, equilibration) must be made
conditional on genuine over/underflow risk. What is a no-op for IEEE can be the dominant
error source for a posit.
