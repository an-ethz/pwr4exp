# Residual Variance-Covariance Matrices

Residual Variance-Covariance Matrices

## Usage

``` r
mkRTrms(data, correlation = NULL)
```

## Arguments

- data:

  a data frame with grouping factors and covariates.

- correlation:

  a `corStruct` object created by
  [corClasses](https://an-ethz.github.io/pwr4exp/reference/corClasses.md)
  functions. If NULL, an identity matrix is assumed.

## Value

A list containing:

- corStruct: An initialized correlation structure.

- corframe: A processed data frame with indexed grouping variables and
  ordering for correlation structures.

- R: A block-diagonal residual variance-covariance structure, not yet
  scaled by the residual variance
