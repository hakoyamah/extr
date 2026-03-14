# Display Modes for Extinction Probabilities

A pair of internal functions for formatting extinction probabilities:
`repr_mode()` decides how to represent a probability \\G\\ based on its
magnitude, and [`format_by_mode()`](format_by_mode.md) formats numbers
accordingly for presentation. These functions affect display only and do
not alter underlying calculations.

## Usage

``` r
repr_mode(g)
```

## Arguments

- g:

  Numeric scalar, a probability in the unit interval. May be non-finite.

## Value

`repr_mode()` returns a character scalar, one of:

- `"linear_G"` – display \\G\\ on the linear scale;

- `"log_G"` – display \\G\\ on log10 scale for tiny probabilities;

- `"log_Q"` – display \\Q=1-G\\ on log10 scale when \\G \approx 1\\;

- `"undefined"` – input was not finite.

[`format_by_mode()`](format_by_mode.md) returns a character scalar
formatted according to the selected mode and the kind of value (point
estimate, lower, or upper CI bound).

## Details

Probabilities very close to 0 or 1 are displayed on a log scale to
improve readability. The thresholds are fixed at `1e-300` ("near zero")
and `1 - 1e-15` ("near one"). Non-finite inputs are reported as
`"undefined"`.

## Author

Hiroshi Hakoyama, <hiroshi.hakoyama@gmail.com>
