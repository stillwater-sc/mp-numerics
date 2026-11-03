# kalman_scalar

Scalar Kalman filter tracking a 1-D position/velocity system across number
types: double, float, posit<32,2>, posit<32,2> with quire accumulation.

## What this measures

- Filter divergence as a function of precision
- Covariance SPD preservation: does P remain positive definite over N steps?
- Innovation consistency: chi-squared test on normalized innovations
- Quire accumulation benefit in the prediction and update dot products

See docs/control.md for the full taxonomy.
