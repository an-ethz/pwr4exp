# Spherical Correlation Structure

Re-exports [nlme::corSpher](https://rdrr.io/pkg/nlme/man/corSpher.html)
from the nlme package.

## Usage

``` r
corSpher(value, form, nugget, metric, fixed)
```

## Arguments

- value:

  numeric value(s) for the parameter(s)

- form:

  a one-sided formula of the form `~ S1 + ... + Sp`, or
  `~ S1 + ... + Sp | g` spatial covariates S1 through Sp and,
  optionally, a grouping factor g.

- nugget:

  logical; if `TRUE` a nugget effect is added

- metric:

  character string specifying the distance metric

- fixed:

  unused

## See also

[corClasses](https://an-ethz.github.io/pwr4exp/reference/corClasses.md)
See [nlme::corSpher](https://rdrr.io/pkg/nlme/man/corSpher.html) for
original documentation.
