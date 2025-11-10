# Building Design Matrices and Covariance Structures for Linear Mixed Models

Constructs design matrices for fixed and random effects, along with
variance-covariance structures for random effects (G-side) and residuals
(R-side).

## Usage

``` r
mkStruct(formula, data, correlation)
```

## Arguments

- formula:

  a model formula.

- data:

  a data frame containing the variables used in the model.

- correlation:

  a `corStruct` object created by
  [corClasses](https://an-ethz.github.io/pwr4exp/reference/corClasses.md)
  functions. If NULL, an identity matrix is assumed.

## Value

A list containing:

- `data`: Processed data frame with NA values omitted.

- `fxTrms`: Fixed-effects design structure, including model frame and
  design matrix.

- `reTrms`: Random-effects structure (if applicable), including grouping
  factors, design matrices, and variance-covariance matrix.

- `rTrms`: Residual structure (R-side variance-covariance components).

- `formula`: Expanded model formula.
