# Control and Stability Numerics

Mixed-precision numerical failures in control algorithms — taxonomy, known
failure modes, and experiment pointers.

## Why control algorithms stress mixed-precision arithmetic

Control algorithms accumulate error differently from batch solvers. A linear
solver runs once; a Kalman filter runs at every timestep, and any rounding
error in the state estimate or covariance matrix compounds over the entire
trajectory. Two properties make this especially acute:

1. **Recursive structure.** Each step's output is the next step's input.
   A representation error that would be negligible in a one-shot solve grows
   without bound if the algorithm is not numerically stable in the Lyapunov sense.

2. **Covariance matrix constraints.** The Kalman filter's error covariance
   matrix P must remain symmetric positive definite (SPD) at every step.
   Finite-precision arithmetic can violate this — P develops negative
   eigenvalues, the filter diverges, and the symptom looks like a control
   failure rather than a numerical one.

## Taxonomy of failure modes

### 1. Representation failures

| Mode | Trigger | Posit behavior |
|------|---------|----------------|
| Overflow | Large state magnitudes or ill-scaled process noise | NaR instead of +-inf; saturates gracefully |
| Underflow | Small covariance entries near convergence | Tapered precision: more bits near zero than IEEE |
| Catastrophic cancellation | Innovation y-Hx when measurement ~= prediction | Quire accumulation eliminates this |
| Denormal drag | IEEE subnormals slow the FPU | Posits have no subnormals; uniform cost |

### 2. Covariance matrix failures

| Mode | Trigger | Detection |
|------|---------|-----------|
| Loss of symmetry | Accumulated rounding in P = FPFt + Q | norm(P - Pt) > eps |
| Loss of positive definiteness | Negative eigenvalue from rounding | Cholesky factorization fails |
| Rank deficiency | Observable/controllable subspace collapses | Condition number of P exceeds 1/eps |

Remedies: Joseph form P = (I-KH)P(I-KH)t + KRKt (symmetric by construction);
UD factorization (Bierman); square-root filter (Carlson/Potter). LDL*
factorization avoids the square-root rounding in Cholesky entirely.

### 3. Accumulation failures in the Kalman update

The innovation covariance S = HPHt + R and the Kalman gain K = PHtS^-1 each
involve dot products of length n. With a 512-bit quire for posit<32,2>, dot
products of length up to 2^31-1 are computed exactly before the single final
rounding -- the "always works" guarantee from The End of Error section 17.4.

### 4. Algorithm instability

| Mode | Description |
|------|-------------|
| Filter divergence | Innovation grows unbounded; gain collapses to zero |
| Consistency loss | Normalized innovation outside expected chi-squared bounds |
| Discretization artifacts | Euler loses stability margins that RK4 or Van Loan preserves |
| Observability rank collapse | States become numerically unestimable |

### 5. Mixed-precision opportunities

The Kalman filter has a natural precision hierarchy:

- Prediction: x- = Fx, P- = FPFt + Q -- quire accumulation beneficial
- Innovation: y~ = y - Hx- -- cancellation risk, quire required
- Gain: K = P-Ht(HP-Ht + R)^-1 -- conditioning of S governs precision needed
- Update: x = x- + Ky~ -- quire beneficial
- Covariance update: P = (I-KH)P- -- Joseph form preferred for SPD

## Experiment pointers

| Experiment | Location | Status |
|-----------|----------|--------|
| Scalar Kalman (1-D position/velocity) | applications/control/kalman_scalar | planned |
| Matrix Kalman precision sweep | benchmarks/control/kalman_precision | planned |
| Covariance SPD preservation | benchmarks/control/spd_preservation | planned |
| Innovation consistency test | benchmarks/control/consistency | planned |

## References

- Gustafson, J.L., The End of Error: Unum Computing, section 17.4
- Bierman, G.J., Factorization Methods for Discrete Sequential Estimation (1977)
- Carlson, N.A., Fast triangular factorization of the square root filter (1973)
- Van Loan, C., Computing integrals involving the matrix exponential (1978)
