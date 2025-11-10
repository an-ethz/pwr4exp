# ARMA(p,q) Correlation Structure

Re-exports [nlme::corARMA](https://rdrr.io/pkg/nlme/man/corARMA.html)
from the nlme package.

## Usage

``` r
corARMA(value, form, p, q, fixed)
```

## Arguments

- value:

  numeric vector of parameter values

- form:

  a one-sided formula of the form `~ t`, `~ 1 | g`, or `~ t | g`,
  specifying a time covariate `t`, a grouping factor `g`, or both. When
  no time covariate is specified, the row order of the data within each
  group is assumed to represent the chronological order of measurements.

- p:

  non-negative integer specifying the autoregressive order

- q:

  non-negative integer specifying the moving average order

- fixed:

  unused

## Details

In the original
[nlme::corARMA](https://rdrr.io/pkg/nlme/man/corARMA.html) function, a
covariate `t` must be an integer class. In pwr4exp, `t` can also be a
factor class, which is then converted to an integer internally for
chronological order. The class of `t` in the model formula, if present,
is not converted. For example, a time covariate can be fitted as a
factor in the model formula, whereas it is converted to an integer in
the correlation formula.

## See also

[corClasses](https://an-ethz.github.io/pwr4exp/reference/corClasses.md)
See [nlme::corARMA](https://rdrr.io/pkg/nlme/man/corARMA.html) for
original documentation.
