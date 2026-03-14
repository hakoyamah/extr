# Confidence Intervals for Extinction Probability (w–z, DI model)

Computes confidence intervals (CIs) for extinction probability in a
density-independent (drifted Wiener) model using the w–z method, and
provides a formatter for display-ready CI strings.

## Usage

``` r
confidence_interval_wz_di(mu, xd, s, th, tq, qq, alpha, prob_fun = ext_prob_di)

ci_wz_format_di(mu, xd, s, th, tq, qq, alpha, digits = 5L)
```

## Arguments

- mu:

  Numeric. Estimated growth rate \\\hat{\mu}\\.

- xd:

  Numeric. Log-distance to threshold \\x_d=\log(n_q/n_e)\\.

- s:

  Numeric. Estimated environmental variance \\\hat{\sigma}^2\\.

- th:

  Numeric. Time horizon \\t^{\ast}\\ for evaluation.

- tq:

  Numeric. Observation span \\t_q\\ (first to last time).

- qq:

  Integer. Number of intervals \\q\\ (sample size minus 1).

- alpha:

  Numeric. Significance level \\\alpha\\\\=\\1\\-\\CL.

- prob_fun:

  Function. One of `ext_prob_di`, `log_ext_prob_di`, `log_ext_comp_di`.

- digits:

  Integer. Significant digits for formatting (used only by
  `ci_wz_format_di`).

## Value

For `confidence_interval_wz_di`: numeric vector `c(lower, upper)` on the
chosen scale (natural log if `log_*`).  
For `ci_wz_format_di`: named character vector `c(lower, upper)` with
values preformatted for display.

## Details

The w–z method derives CIs by inverting noncentral-\\t\\ distributions
for the transformed statistics \\w\\ and \\z\\, then combining them to
form bounds on \\G(w,z)\\.

Exact confidence intervals for \\w\\ and \\z\\ are obtained by
numerically solving the noncentral-\\t\\ quantile equations
corresponding to the observed statistics.

The CI for \\G(w,z)\\ is then approximated by \$\$ \bigl(
G(\overline{w},\\\underline{z}),\\ G(\underline{w},\\\overline{z})
\bigr), \$\$ where \\\overline{w}\\, \\\underline{w}\\,
\\\overline{z}\\, and \\\underline{z}\\ are the upper and lower
confidence limits for \\w\\ and \\z\\, respectively. Across the full
parameter space, this approach achieves near-nominal coverage. When
\$z\$ is large and positive, \$G\$ depends only on \$w\$, so exact CIs
are available.

The argument `prob_fun` selects which probability is evaluated:

- `ext_prob_di`, `log_ext_prob_di`: CIs for \\G(w,z)\\ and \\\log
  G(w,z)\\; returned as (lower, upper).

- `log_ext_comp_di`: CIs for \\\log Q(w,z)= \log (1-G(w,z))\\; returned
  as (lower, upper).

## Author

Hiroshi Hakoyama, <hiroshi.hakoyama@gmail.com>
