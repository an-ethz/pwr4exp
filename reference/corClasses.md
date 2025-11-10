# Correlation Structure Classes

Standard classes of correlation structures (`corStruct`) available from
the nlme package and re-exported by pwr4exp for convenience when
specifying correlation structures in
[`mkdesign`](https://an-ethz.github.io/pwr4exp/reference/mkdesign.md).

All arguments are identical to the corresponding nlme functions. For
more details on the original implementations, see
[nlme::corClasses](https://rdrr.io/pkg/nlme/man/corClasses.html).

Note: In the original
[nlme::corAR1](https://rdrr.io/pkg/nlme/man/corAR1.html),
[nlme::corARMA](https://rdrr.io/pkg/nlme/man/corARMA.html), and
[nlme::corSymm](https://rdrr.io/pkg/nlme/man/corSymm.html) functions,
the covariate `t` in the correlation formula `~ t` or `~ t | g` must be
an integer class. In pwr4exp, the covariate can also be a factor class,
which is then converted to an integer internally for sorting purposes.
The class of the covariate variable in the model formula, if present,
will not be converted. For example, a time covariate can be fitted as a
factor in the model formula, whereas it is converted to an integer in
the correlation formula temporarily for matrix sorting.

## Details

Available standard classes:

- [`corAR1`](https://an-ethz.github.io/pwr4exp/reference/corAR1.md):

  autoregressive process of order 1.

- [`corARMA`](https://an-ethz.github.io/pwr4exp/reference/corARMA.md):

  autoregressive moving average process, with arbitrary orders for the
  autoregressive and moving average components.

- [`corCAR1`](https://an-ethz.github.io/pwr4exp/reference/corCAR1.md):

  continuous AR(1)

- [`corCompSymm`](https://an-ethz.github.io/pwr4exp/reference/corCompSymm.md):

  compound symmetry structure corresponding to a constant correlation.

- [`corExp`](https://an-ethz.github.io/pwr4exp/reference/corExp.md):

  exponential spatial correlation.

- [`corGaus`](https://an-ethz.github.io/pwr4exp/reference/corGaus.md):

  Gaussian spatial correlation.

- [`corLin`](https://an-ethz.github.io/pwr4exp/reference/corLin.md):

  linear spatial correlation.

- [`corRatio`](https://an-ethz.github.io/pwr4exp/reference/corRatio.md):

  Rational quadratics spatial correlation.

- [`corSpher`](https://an-ethz.github.io/pwr4exp/reference/corSpher.md):

  spherical spatial correlation.

- [`corSymm`](https://an-ethz.github.io/pwr4exp/reference/corSymm.md):

  general correlation matrix, with no additional structure.

## References

Pinheiro, J.C., and Bates, D.M. (2000) "Mixed-Effects Models in S and
S-PLUS", Springer.
