# Computes the w and z Statistics

Estimators for the \\w\\ and \\z\\ statistics used in extinction
probability calculations under a density-independent population model.

## Usage

``` r
w_statistic(mu, xd, s, th)

z_statistic(mu, xd, s, th)
```

## Arguments

- mu:

  numeric: Estimated population growth rate, \\\hat{\mu}\\.

- xd:

  numeric: Distance to extinction threshold on a log scale, \\x_d =
  \log(n_q / n_e)\\.

- s:

  numeric: Estimated environmental variance, \\\hat{\sigma}^2\\.

- th:

  numeric: Time horizon for extinction probability evaluation, denoted
  \\t^{\ast}\\.

## Value

numeric: Value of the statistic.

## Details

The statistics are defined as \$\$ \hat w = \frac{\hat \mu t^{\ast} +
x_d}{\sqrt{\hat \sigma^2 t^{\ast}}}, \qquad \hat z = \frac{- \hat \mu
t^{\ast} + x_d}{\sqrt{\hat \sigma^2 t^{\ast}}}. \$\$

## Author

Hiroshi Hakoyama, <hiroshi.hakoyama@gmail.com>
