# Version 1.1.0

This is an update to a previously submitted CRAN package.

## Changes in this version

* Added `method = "oear"` to `ext_di()`, allowing extinction probability estimates and associated confidence intervals to be computed using an observation-error-and-autocovariance-robust (OEAR) plug-in variance estimate.
* Implemented the OEAR variance estimator using HAC long-run variance estimation with AR(1) pre-whitening and a Bartlett kernel.
* Added OEAR-related diagnostics to the returned object and printed output.
* Updated documentation, tests, and package metadata accordingly.

## Test environments

* macOS (Monterey 12.7.6), R 4.5.2, aarch64-apple-darwin20 — `R CMD check --no-manual --as-cran` clean
* macOS (Sonoma 14.8.2), R 4.5.2, aarch64-apple-darwin20 — `R CMD check --no-manual --as-cran` clean
* Linux (Clear Linux OS), R 4.5.2, x86_64-pc-linux-gnu — `R CMD check --no-manual` clean
* Windows Server 2022 x64, R 4.5.2 Patched (2026-02-13 r89426 ucrt), x86_64-w64-mingw32 — win-builder release check clean
* Windows Server 2022 x64, R Under development (unstable) (2026-03-08 r89578 ucrt), x86_64-w64-mingw32 — win-builder devel check clean

## Reverse dependencies

There are currently no reverse dependencies.
