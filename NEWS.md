# extr 1.0.0 (2025-09-15)

## Up to 0.9.0
- Core extinction-risk estimator `ext_di` for density-independent
  (drifted Wiener) models. Provides:
  - MLEs of growth rate and environmental variance from irregular time series.
  - Point estimate of extinction probability.
  - Confidence intervals for extinction probability (w-z method) with
    near-nominal coverage.
  - QQ-plot diagnostic optional.

## Changes in 1.0.0
- Numerical stability and range
  - Adaptive probability representation on G, log G, and log(1-G). The
    complement is used when G is near 0 or 1 to maintain accuracy.
  - The log scale removes double-precision bounds near 0 and 1, extending
    the safe window down to symbolic exp(-DBL_MAX). Probabilities close to
    0 or 1 are displayed as 10^-* or 1-10^-*.

- Printed output
  - Printed precision is configurable via options("extr.digits").
