# General Correlation Structure

Re-exports [nlme::corSymm](https://rdrr.io/pkg/nlme/man/corSymm.html)
from the nlme package.

## Usage

``` r
corSymm(value, form, fixed)
```

## Arguments

- value:

  numeric vector of correlation parameter values

- form:

  a one-sided formula of the form `~ t`, `~ 1 | g`, or `~ t | g`,
  specifying a covariate `t`, a grouping factor `g`, or both. A
  covariate indicates the order of the data rows in the residual error
  matrix of each group.

- fixed:

  unused

## See also

[corClasses](https://an-ethz.github.io/pwr4exp/reference/corClasses.md)
See [nlme::corSymm](https://rdrr.io/pkg/nlme/man/corSymm.html) for
original documentation.
