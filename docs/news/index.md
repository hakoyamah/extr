# Changelog

## extr 1.1.0 (2026-03-14)

### Changes in 1.1.0

- Added `method = "oear"` to [`ext_di()`](../reference/ext_di.md),
  allowing extinction probability estimates and associated confidence
  intervals to be computed using an
  observation-error-and-autocovariance-robust (OEAR) plug-in variance
  estimate.
- Implemented the OEAR variance estimator using HAC long-run variance
  estimation with AR(1) pre-whitening and a Bartlett kernel.
- Added OEAR diagnostics (`rho_tilde_pw` and `j`) to the returned object
  and printed output.
- Updated method-specific output: for `method = "naive"`, output
  includes `Variance` and `aic`; for `method = "oear"`, output includes
  `Process variance (OEAR)`, reports `aic` as `NA`, and reports
  `lower_cl_s` and `upper_cl_s` using a chi-square-based pragmatic
  plug-in approximation.

## extr 1.0.0 (2025-09-15)

CRAN release: 2025-09-21

### Up to 0.9.0

- Core extinction-risk estimator `ext_di` for density-independent
  (drifted Wiener) models.
- Includes MLEs of growth rate and environmental variance from irregular
  time series.
- Includes point estimate of extinction probability.
- Includes confidence intervals for extinction probability (w-z method)
  with near-nominal coverage.
- Includes optional QQ-plot diagnostic.

### Changes in 1.0.0

- Numerical stability and range:
  - Adaptive probability representation on G, log G, and log(1-G).
  - Uses complement form when G is near 0 or 1 to maintain numerical
    accuracy.
  - Log-scale representation extends safe range down to symbolic
    exp(-DBL_MAX).
  - Probabilities near 0 or 1 are displayed as `10^-*` or `1-10^-*`.
- Printed output:
  - Printed precision is configurable via `options("extr.digits")`.
