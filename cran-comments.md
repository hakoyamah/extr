# Version 1.0.0

Initial CRAN submission.

## Test environments

* Local
  - macOS (Sonoma 14.7.7), R 4.5.1, aarch64-apple-darwin20 — R CMD check clean
  - Linux (Ubuntu 22.04.5 LTS), R 4.5.1, x86_64-pc-linux-gnu — R CMD check clean

* Remote
  - Windows (release): devtools::check_win_release()
    - Windows Server 2022 x64, R 4.5.1, x86_64-w64-mingw32 — clean
  - Windows (devel): devtools::check_win_devel()
    - Windows Server 2022 x64, R Under development (unstable)
      (2025-09-12 r88822 ucrt), x86_64-w64-mingw32 — clean

## devtools::check() (local, source)
0 errors | 0 warnings | 0 notes

## R CMD check --as-cran (local, source tarball)
Status: 1 NOTE

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Hiroshi Hakoyama <hiroshi.hakoyama@gmail.com>'
  New submission

This NOTE is expected for an initial submission and is not specific to
the package.

## Reverse dependencies

There are currently no reverse dependencies.
