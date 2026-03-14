# Print Extinction Risk Estimates

S3 method that prints a formatted summary for an `ext_di` object. The
output includes the plug-in extinction probability and its CI,
method-specific parameter estimates, a brief data summary, and input
settings.

## Usage

``` r
# S3 method for class 'ext_di'
print(x, digits = NULL, ...)
```

## Arguments

- x:

  An object of class `"ext_di"` as returned by [`ext_di`](ext_di.md).

- digits:

  Integer. Significant digits to display; defaults to `x$digits` when
  `NULL`.

- ...:

  Additional arguments (currently ignored).

## Value

Invisibly returns `x` after printing the formatted results.
