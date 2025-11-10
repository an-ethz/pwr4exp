# AR(1) Correlation Structure

Re-exports [nlme::corAR1](https://rdrr.io/pkg/nlme/man/corAR1.html) from
the nlme package.

## Usage

``` r
corAR1(value, form, fixed)
```

## Arguments

- value:

  numeric value for the correlation parameter.

- form:

  a one-sided formula of the form `~ t`, `~ 1 | g`, or `~ t | g`,
  specifying a time covariate `t`, a grouping factor `g`, or both. When
  no time covariate is specified, the row order of the data within each
  group is assumed to represent the chronological order of measurements.

- fixed:

  unused

## Details

In the original [nlme::corAR1](https://rdrr.io/pkg/nlme/man/corAR1.html)
function, a covariate `t` must be an integer class. In pwr4exp, `t` can
also be a factor class, which is then converted to an integer internally
for chronological order. The class of `t` in the model formula, if
present, is not converted. For example, a time covariate can be fitted
as a factor in the model formula, whereas it is converted to an integer
in the correlation formula.

## See also

[corClasses](https://an-ethz.github.io/pwr4exp/reference/corClasses.md)
See [nlme::corAR1](https://rdrr.io/pkg/nlme/man/corAR1.html) for
original documentation.
