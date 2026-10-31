# Givens rotations and Householder reflectors, one at a time

mp-numerics#13, part of the series in stillwater-sc/mtl5#666. This study measures
how a single plane rotation or a single Householder reflector behaves in each
number system, before any factorization is built on them. It tests MTL5's own
kernels, `mtl::apply_givens_rotation` (`operation/givens.hpp`) and
`mtl::householder` (`operation/householder.hpp`), against the forms LAPACK uses
and against two power-of-two-scaled forms written here.

This is the first of two parts. It covers the hardware IEEE types, posits and
fixed point. The second part adds the cfloat matrix (fp8 to fp256, IEEE-like,
without subnormals, and saturating), lns, takum and dd/qd, together with plots
and the findings note for the MTL5 docs.

## Running it

```bash
cmake --build build --target givens_householder_study
./build/benchmarks/linalg/givens_householder/givens_householder_study results/
python3 benchmarks/linalg/givens_householder/summarize.py results/
```

`--quick` runs a small sample, and `--family ieee,posit,fixpnt` picks families.
The full run takes under a minute. Each CSV has a `.sysinfo` sidecar recording
the platform, compiler, sample counts and the commits of mp-numerics, MTL5 and
Universal. The commits are read when CMake configures, so reconfigure before a
run whose results are kept; `+dirty` marks uncommitted changes. The inputs come
from `std::mt19937_64` and `std::cos`/`std::sin`, so another platform's libm can
move an input by an ulp and change the last digits.

## What is measured

Each operator is built for x = (a, b) and measured exactly. The computed values
are read out of the number type as exact Shewchuk expansions
(`include/sw/mp_numerics/exact.hpp`), so every error below is the error of the
operator the arithmetic produced, rounded once to double for the CSV. Values
beyond double's range are first scaled by a power of two inside the type.

| column | meaning |
|---|---|
| `orth` | ‖QᵀQ − I‖₂ of the computed operator: \|c² + s² − 1\| for a rotation, \|g\|·\|g − 2\| with g = τ vᵀv for a reflector |
| `r_err` | \|stored entry − ‖x‖\| / ‖x‖, for the value a factorization keeps in place of a (r, or the first entry of H x) |
| `annih` | \|second entry of Q x\| / ‖x‖, with Q applied exactly |
| `annih_fp` | the same entry computed in the type |
| `nonfinite` | samples with an infinity or NaN in the operator or its result |
| `identity` | samples where the routine returned the identity |
| `gross` | samples with `orth` or `annih` above 0.5: the operator is wrong. `r_err` is left out, because at the ends of a tapered or fixed-point range it mostly measures how coarsely the type represents ‖x‖ |

`eps` is `std::numeric_limits<T>::epsilon()`, at 1 for a posit; the tables
below give errors in units of it.

### The operators

| name | what it is |
|---|---|
| `givens_mtl5` | `mtl::apply_givens_rotation`: r = sqrt(a·a + b·b), c = a/r, s = b/r. GMRES uses it, and MINRES and the symmetric QR iteration use the same formula |
| `givens_max_scaled` | divide by m = max(\|a\|, \|b\|) first, as LAPACK's `lapy2` does; c = a/r |
| `givens_safe_scaled` | Anderson's safe scaling, transcribed from Reference-LAPACK `dlartg.f90` (3.10), its thresholds taken from `std::numeric_limits` |
| `givens_pow2_scaled` | scale by an exact 2^-e with 2^e ≤ max(\|a\|, \|b\|), so the squares sit near 1; c and s from the scaled values |
| `householder_mtl5` | `mtl::householder`: max-scaled, Parlett's sign (onto +‖x‖e₁ for x₀ > 0), v(0) = 1, applied as `qr_solve` applies it |
| `householder_textbook` | `mtl::householder(x, tag::textbook_reflector)`: v = x − αe₁ with α = −sign(x₀)‖x‖, unscaled, unnormalized |
| `householder_lapack_pow2` | `dlarfg`'s sign (onto −sign(x₀)‖x‖e₁) and v(0) = 1, on x scaled by an exact 2^-e |

### The inputs

| input | x | per row |
|---|---|---|
| `base` | random directions, ‖x‖ = 1 | 10000 samples |
| `scale` | random directions, ‖x‖ = 2^k for k across the type's whole range (`param` = k) | 32 samples |
| `tail` | (±1, ±2^-j (1 + u)), j from 0 to 8 past the type's precision (`param` = j) | 32 samples |
| `repeat` | `power`: one operator applied k times to 16 unit vectors (`dev` = max \| ‖y‖ − 1 \|); `product`: k fresh operators accumulated into an 8 × 8 Q (`dev` = ‖QᵀQ − I‖_F) | 8 and 4 trials |

The types: `float`, `double`, `long double` (x87 80-bit here; the `bits`
column records which one ran), posit<8,2>, <16,2>, <32,2>, <64,2> and
posit<16,1>, and fixpnt<16,8> and <32,16>, each wrapping (modulo) and
saturating.

## Findings

**1. MTL5's rotation fails outside half the exponent range.** a·a overflows
once \|a\| passes sqrt(max), and underflows below sqrt(min), long before r
would. The clean band is k ∈ [−74, 63] for float out of [−149, 127], and
[−533, 509] for double. How it fails depends on the type. In IEEE, r becomes
infinite, so c = s = 0 and the rotation is the zero matrix; or r underflows to
0, and the rotation is the identity with r = 0 stored for ‖x‖. In a posit or a
saturating fixed-point type, the square saturates and the rotation is finite
and wrong. A wrapping fixed-point type wraps. Even inside the band, a posit's
squares land where it has few bits: posit<64,2> reaches 5×10⁸ eps.

**2. `lapy2`'s scaling stops the overflow, but not the loss of accuracy.**
With c = a/r, the error of r at the input's own magnitude goes into c. That
magnitude is subnormal in IEEE, or far from 1 in a posit, where r has few bits.
The operator is not grossly wrong there, but its orthogonality error reaches
0.25: 1.1×10¹⁵ eps in double near the subnormal range, and 1.4×10¹⁷ eps in
posit<64,2>.

**3. `dlartg`'s safe scaling guards the range, not the precision.** In IEEE
it is what it was designed to be: clean across the whole range, within 2.1 eps.
In a posit it never fails grossly, but its thresholds come from the range: it
scales only when an entry is beyond sqrt(min) or sqrt(max/2), and inside those
bounds a posit's squares fall where it has few bits. posit<64,2> reaches
5×10⁸ eps and posit<32,2> 9×10³ eps. The thresholds would have to come from
where the type is precise, not from where it overflows.

**4. Scaling by an exact power of two fixes the rotation in every type.**
`givens_pow2_scaled` is clean across the whole range of every type, within
2.2 eps, posits and fixed point included. The scaling rounds nothing in a
floating-point type and puts the squares near 1.

**5. MTL5's reflector loses accuracy on small tails in tapered and
fixed-point types, because of its choice of sign.** Parlett's sign gives
v₀ = −σ/(x₀ + ‖x‖), which is small when the tail is small, so after
normalizing v(0) = 1 the second entry is v₁ ≈ −2x₀/x₁, which is large. In IEEE
that costs nothing (about 9 eps at worst). In a posit, v₁ and 1 + v₁² land where the
posit has few bits: on the tail sweep posit<32,2> reaches 512 eps and
posit<64,2> 8×10⁵ eps. In fixed point, which has absolute precision, the
reflector reaches 2.4×10³ eps in fixpnt<32,16>. In the tail sweep x₀ = ±1,
so MTL5's scaling divides by 1 and rounds nothing: the sign alone is the
cause. With `dlarfg`'s sign, |v₀| = |x₀| + ‖x‖ ≥ ‖x‖, so every entry of the
normalized v is at most 1, and `householder_lapack_pow2` stays within 5.4 eps
for every type, at every scale and every tail. MTL5's identity fallback
(mtl5#662) catches the worst cases, but only once the error exceeds the tail:
from j = 8 in posit<16,2> and j = 40 in posit<64,2>.

**6. MTL5's textbook reflector is accurate but unscaled.** It has the LAPACK
sign and no normalization, so it is within 6.4 eps on every tail, and within
1 eps in fixed point. Its unscaled sum of squares limits it to the same
half-range band as the unscaled rotation.

**7. Repeated application.** One rounded rotation applied k times shrinks or
stretches every vector by about k·\|c² + s² − 1\|/2, so the error grows
linearly: 1.4×10³ to 4.7×10³ eps at k = 4096 in every type with enough
precision, whichever way the rotation was computed. Accumulating k fresh
operators into an 8 × 8 Q grows far more slowly: 40 to 250 eps at k = 1024 in
IEEE and in posits from 16 bits up, and up to 510 eps in fixed point.
posit<8,2> is too narrow for either.

### What this suggests for MTL5

These are candidates, not changes made here:

- `apply_givens_rotation`, and the same formula in MINRES and the symmetric
  QR iteration, could scale by an exact power of two before squaring, as
  `givens_pow2_scaled` does. In IEEE that costs two exponent extractions and
  removes the half-range limit; in a posit or fixed point it also keeps the
  accuracy.
- `householder()` could reflect onto −sign(x₀)‖x‖e₁, as `dlarfg` does. That
  changes the sign of R's diagonal, which MTL5 currently keeps positive.

## Tables

Generated by `summarize.py` from `results/`. Each entry is in units of the
type's eps.

### Scale sweep

Each entry gives the range of k with no gross or non-finite sample, then the
worst orthogonality within it.

| type | range of k | `givens_mtl5` | `givens_max_scaled` | `givens_safe_scaled` | `givens_pow2_scaled` | `householder_mtl5` | `householder_textbook` | `householder_lapack_pow2` |
|---|---|---|---|---|---|---|---|---|
| float | -149..127 | -74..63, 3.89 | -147..127, 2.1e+06 | all, 1.99 | all, 1.99 | all, 7.62 | -74..63, 3.5e+06 | all, 4.42 |
| double | -1074..1023 | -533..509, 1.23 | -1072..1023, 1.1e+15 | all, 1.85 | all, 1.46 | all, 8.44 | -533..509, 1.9e+13 | all, 4.31 |
| long double | -16445..16383 | -8186..8124, 1.99 | -16443..16383, 2.3e+18 | all, 1.99 | all, 2.07 | all, 9.09 | -8186..8124, 3.22 | all, 4.76 |
| posit<8,2> | -24..24 | -12..12, 4 | -11..17, 2 | -24..24, 4 | all, 1.32 | all, 3.47 | -4..1, 2.31 | all, 3.22 |
| posit<16,2> | -56..56 | -28..28, 72 | -43..44, 512 | all, 42.5 | all, 1.8 | all, 15.6 | -18..17, 682 | all, 3.83 |
| posit<32,2> | -120..120 | -60..60, 2.0e+04 | -108..108, 3.4e+07 | all, 9.08e+03 | all, 1.98 | all, 69.6 | -49..49, 4.4e+07 | all, 4.76 |
| posit<64,2> | -248..248 | -120..120, 5.0e+08 | -232..239, 1.4e+17 | all, 5.0e+08 | all, 2.01 | all, 33.6 | -114..111, 1.9e+17 | all, 4.28 |
| posit<16,1> | -28..28 | -13..14, 45 | -25..24, 1.02e+03 | all, 45 | all, 1.44 | all, 68.6 | -10..9, 1.36e+03 | all, 3.76 |
| fixpnt<16,8,mod> | -8..6 | -4..3, 16 | -7..6, 64 | all, 12.2 | all, 1.63 | all, 67.4 | -3..2, 59.4 | all, 2.94 |
| fixpnt<16,8,sat> | -8..6 | -4..3, 50.6 | -7..6, 64 | all, 12.2 | all, 1.56 | all, 67.4 | -3..2, 59.4 | all, 2.75 |
| fixpnt<32,16,mod> | -16..14 | -8..7, 327 | -15..14, 1.6e+04 | all, 252 | all, 1.83 | all, 2.01e+03 | -7..6, 1.3e+04 | all, 4.43 |
| fixpnt<32,16,sat> | -16..14 | -8..7, 1.03e+03 | -15..14, 1.6e+04 | all, 1.03e+03 | all, 1.83 | all, 2.08e+03 | -7..6, 2.5e+04 | all, 3.68 |

### Tail sweep

The worst orthogonality, and the first j at which the routine returns the
identity (I).

| type | `givens_mtl5` | `givens_max_scaled` | `givens_safe_scaled` | `givens_pow2_scaled` | `householder_mtl5` | `householder_textbook` | `householder_lapack_pow2` |
|---|---|---|---|---|---|---|---|
| float | 1.61 | 1.61 | 1.61 | 1.61 | 8.3 (I from j=22) | 3.49 | 4.07 |
| double | 1.77 | 1.66 | 1.77 | 1.77 | 7.43 (I from j=51) | 3.89 | 3.98 |
| long double | 1.66 | 2.26 | 1.66 | 1.66 | 8.53 (I from j=62) | 3.85 | 3.71 |
| posit<8,2> | 1.32 | 1.69 | 1.32 | 1.32 | 3.47 (I from j=3) | 2.53 | 3.22 |
| posit<16,2> | 1.88 | 1.48 | 1.88 | 1.88 | 23.3 (I from j=8) | 3.66 | 4.02 |
| posit<32,2> | 1.79 | 2.22 | 1.79 | 1.79 | 512 (I from j=18) | 3.56 | 3.8 |
| posit<64,2> | 1.91 | 1.91 | 1.91 | 1.91 | 8.0e+05 (I from j=40) | 3.9 | 4.14 |
| posit<16,1> | 1.59 | 2.08 | 1.59 | 1.59 | 56.1 (I from j=6) | 6.31 | 3.63 |
| fixpnt<16,8,mod> | 1.44 | 2.15 | 1.44 | 1.44 | 63.8 (I from j=3) | 0.738 | 2.3 (I from j=10) |
| fixpnt<16,8,sat> | 1.67 | 1.86 | 1.67 | 1.67 | 67.4 (I from j=2) | 0.795 | 2.86 (I from j=10) |
| fixpnt<32,16,mod> | 1.77 | 1.77 | 1.77 | 1.77 | 2.36e+03 (I from j=5) | 0.873 | 2.74 (I from j=18) |
| fixpnt<32,16,sat> | 1.99 | 2.6 | 1.99 | 1.99 | 2.41e+03 (I from j=5) | 0.926 | 2.78 (I from j=18) |

The fixed-point identities of `householder_lapack_pow2` start where the tail
rounds to zero in the type.
