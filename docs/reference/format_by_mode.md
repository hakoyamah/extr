# Format a Probability using a Display Mode

Format a probability (point estimate or CI bound) according to a chosen
representation mode. Allowed modes are `"linear_G"`, `"log_G"`,
`"log_Q"`, and `"undefined"`. The function performs no inference; it
only formats the numeric values already computed (linear and log).

Notes on arguments:

- `log_g`, `log_q`, `ci_log_g`, `ci_log_q` are natural logarithms (from
  [`log()`](https://rdrr.io/r/base/Log.html) in R).

- For `"log_Q"`: the upper bound for \\G\\ uses `ci_log_q[1]` and the
  lower bound uses `ci_log_q[2]`. This ordering reflects that \\Q =
  1-G\\ decreases as \\G\\ increases.

## Usage

``` r
format_by_mode(
  mode,
  kind,
  linear_g,
  log_g,
  log_q,
  ci_linear_g,
  ci_log_g,
  ci_log_q,
  digits
)
```

## Arguments

- mode:

  Character scalar. One of `"linear_G"`, `"log_G"`, `"log_Q"`,
  `"undefined"`.

- kind:

  Character scalar. One of `"point"`, `"lower"`, `"upper"`.

- linear_g:

  Numeric scalar. Point estimate on the linear scale.

- log_g:

  Numeric scalar. `log(G)` (natural log). May be non-finite.

- log_q:

  Numeric scalar. `log(Q)` with \\Q=1-G\\ (natural log). May be
  non-finite.

- ci_linear_g:

  Numeric length-2: `c(lower, upper)` for \\G\\.

- ci_log_g:

  Numeric length-2: `c(lower_logG, upper_logG)`.

- ci_log_q:

  Numeric length-2: `c(upper_logQ, lower_logQ)`.

- digits:

  Integer. Significant digits for numeric formatting.

## Value

A character scalar formatted for display. Returns `"NA"` if the required
value is non-finite or missing.

## Author

Hiroshi Hakoyama, <hiroshi.hakoyama@gmail.com>
