# Computes Extinction Probability for a Density-Independent Model

Compute \\G(w,z)\\ and its complement \\Q(w,z)=1-G(w,z)\\ for a
density-independent (drifted Wiener) model, on both linear and log
scales.

## Usage

``` r
ext_prob_di(w, z)

log_ext_prob_di(w, z)

log_ext_comp_di(w, z)

ext_prob_format_di(w, z, digits = 5L)
```

## Arguments

- w:

  Numeric; transformed parameter \\w=(\mu t+x_d)/(\sigma\sqrt{t})\\.

- z:

  Numeric; transformed parameter \\z=(-\mu t+x_d)/(\sigma\sqrt{t})\\.

- digits:

  Integer; significant digits for formatting (only for
  `ext_prob_format_di`).

## Value

For `ext_prob_di`: numeric scalar \\G(w,z)\\.  
For `log_ext_prob_di`: numeric scalar \\\log G(w,z)\\.  
For `log_ext_comp_di`: numeric scalar \\\log Q(w,z)\\.  
For `ext_prob_format_di`: character string formatted for display.

## Details

For any \\t\>0\\ with \\w+z\>0\\, \$\$ \Pr\[T \leq t\] =
G(w,z)=\Phi(-w)+ \exp\\\left(\tfrac{z^2-w^2}{2}\right)\Phi(-z), \qquad
Q(w,z)=1-G(w,z). \$\$ Here \\\Phi\\ and \\\phi\\ denote the standard
normal CDF and PDF.

**Stability strategy.** (i) For large \\z\\, rewrite the product
\\\exp((z^2-w^2)/2)\\\Phi(-z)\\ via the Mills ratio and replace it by an
8-term asymptotic series when \\z \ge 19\\. (ii) On the log scale, use
log-sum-exp and a stable log-difference (`log1mexp`) built from
`log1p`/`expm1` to retain tail info.

**Domain.** Scalar inputs are assumed and require \\w+z\>0\\.

**Functions.**

- `ext_prob_di(w,z)`: returns \\G(w,z)\\ (linear scale).

- `log_ext_prob_di(w,z)`: returns \\\log G(w,z)\\.

- `log_ext_comp_di(w,z)`: returns \\\log Q(w,z)\\.

- `ext_prob_format_di(w,z,digits)`: formats a point estimate using
  [`repr_mode()`](repr_mode.md) and
  [`format_by_mode()`](format_by_mode.md).

## See also

[`statistics_di`](statistics_di.md), [`repr_mode`](repr_mode.md),
[`format_by_mode`](format_by_mode.md)

## Author

Hiroshi Hakoyama, <hiroshi.hakoyama@gmail.com>
