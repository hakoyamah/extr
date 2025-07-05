
<!-- README.md is generated from README.Rmd. Please edit that file -->

# extr

The goal of **extr** is to provide tools for estimating extinction
probabilities and confidence intervals for stochastic population models.

## Installation

``` r
install.packages("extr_1.0.0.tar.gz", repos = NULL, type="source")

# install.packages("remotes") # if needed
remotes::install_github("hakoyamah/extr")
```

## Function

### `ext_di()`

Estimates demographic parameters and extinction probabilities from time
series of population sizes based on a Wiener-drift process model.
Implements the w-z method (Hakoyama, 2025) to improve the accuracy of
confidence intervals over conventional approaches

## Example

``` r
library(extr)
```

``` r
dat <- data.frame(Time = c(1959, 1960, 1961, 1962, 1963, 1964, 1965,
1966, 1967, 1968, 1969, 1970, 1971, 1972, 1973, 1974, 1975, 1976,
1977, 1978, 1979, 1980, 1981, 1982, 1983, 1984, 1985, 1986, 1987),
Population = c(44, 47, 46, 44, 46, 45, 46, 40, 39, 39, 42, 44, 41, 40,
33, 36, 34, 39, 35, 34, 38, 36, 37, 41, 39, 51, 47, 57, 47))
```

``` r
ext_di(dat, th = 100)
#> --- Estimates ---
#>                                                       Estimate
#> Probability of decline to 1 within 100 years (MLE): 9.4128e-05
#> Growth rate (MLE):                                   0.0023556
#> Variance (MLE):                                        0.01087
#> Unbiased variance:                                    0.011273
#> AIC for the distribution of N:                          165.06
#>                                                                        CI
#> Probability of decline to 1 within 100 years (MLE):  (1.4586e-13, 0.5653)
#> Growth rate (MLE):                                  (-0.038814, 0.043525)
#> Variance (MLE):                                     (0.0070464, 0.020885)
#> Unbiased variance:                                                      -
#> AIC for the distribution of N:                                          -
#> 
#> --- Data Summary ---
#>                               Value
#> Current population size, nq:     47
#> xd = ln(nq / ne):            3.8501
#> Sample size, q + 1:              29
#> 
#> --- Input Parameters ---
#>                                                         Parameter
#> Time unit:                                                  years
#> Extinction threshold of population size, ne:                    1
#> Time window for extinction risk evaluation (years), th:       100
#> Significance level, alpha:                                   0.05
```

``` r
ext_di(dat, th = 100, ne = 10)
#> --- Estimates ---
#>                                                       Estimate
#> Probability of decline to 10 within 100 years (MLE):  0.096852
#> Growth rate (MLE):                                   0.0023556
#> Variance (MLE):                                        0.01087
#> Unbiased variance:                                    0.011273
#> AIC for the distribution of N:                          165.06
#>                                                                         CI
#> Probability of decline to 10 within 100 years (MLE):  (1.0698e-05, 0.9898)
#> Growth rate (MLE):                                   (-0.038814, 0.043525)
#> Variance (MLE):                                      (0.0070464, 0.020885)
#> Unbiased variance:                                                       -
#> AIC for the distribution of N:                                           -
#> 
#> --- Data Summary ---
#>                               Value
#> Current population size, nq:     47
#> xd = ln(nq / ne):            1.5476
#> Sample size, q + 1:              29
#> 
#> --- Input Parameters ---
#>                                                         Parameter
#> Time unit:                                                  years
#> Extinction threshold of population size, ne:                   10
#> Time window for extinction risk evaluation (years), th:       100
#> Significance level, alpha:                                   0.05
```

``` r
ext_di(dat, th = 100, ne = 10, qq_plot = TRUE)
```

<img src="man/figures/README-qqplot-1.png" width="100%" />

    #> --- Estimates ---
    #>                                                       Estimate
    #> Probability of decline to 10 within 100 years (MLE):  0.096852
    #> Growth rate (MLE):                                   0.0023556
    #> Variance (MLE):                                        0.01087
    #> Unbiased variance:                                    0.011273
    #> AIC for the distribution of N:                          165.06
    #>                                                                         CI
    #> Probability of decline to 10 within 100 years (MLE):  (1.0698e-05, 0.9898)
    #> Growth rate (MLE):                                   (-0.038814, 0.043525)
    #> Variance (MLE):                                      (0.0070464, 0.020885)
    #> Unbiased variance:                                                       -
    #> AIC for the distribution of N:                                           -
    #> 
    #> --- Data Summary ---
    #>                               Value
    #> Current population size, nq:     47
    #> xd = ln(nq / ne):            1.5476
    #> Sample size, q + 1:              29
    #> 
    #> --- Input Parameters ---
    #>                                                         Parameter
    #> Time unit:                                                  years
    #> Extinction threshold of population size, ne:                   10
    #> Time window for extinction risk evaluation (years), th:       100
    #> Significance level, alpha:                                   0.05
