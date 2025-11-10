# Design Matrices and Variance Components for Random Effects

Adapted from lme4, this function constructs the design matrix (Z),
variance-covariance matrix (G), etc.

## Usage

``` r
mkReTrms(
  bars,
  fr,
  drop.unused.levels = TRUE,
  reorder.terms = FALSE,
  reorder.vars = FALSE
)
```

## Arguments

- bars:

  a list of parsed random-effects terms

- fr:

  a model frame in which to evaluate these terms

- drop.unused.levels:

  (logical) drop unused factor levels?

- reorder.terms:

  arrange random effects terms in decreasing order of number of groups
  (factor levels)?

- reorder.vars:

  arrange columns of individual random effects terms in alphabetical
  order?

## Value

A list with the following components:

- Zt: Transposed random-effects design matrix.

- G: Variance-covariance matrix for random effects.

- Gind: Index mapping variance-covariance parameters to their positions
  in G.

- G_temp: List of individual variance component matrices for each random
  effect.

- flist: List of grouping factors used in random effects.

- cnms: Column names of the random-effects design matrix.

- Ztlist: List of per-term transposed design matrices for random
  effects.

- nl: Number of levels for each grouping factor.
