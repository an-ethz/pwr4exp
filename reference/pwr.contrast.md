# Power of contrasts

Computes the statistical power of t-tests for comparisons among means.

## Usage

``` r
pwr.contrast(
  object,
  which,
  by = NULL,
  contrast = c("pairwise", "poly", "trt.vs.ctrl"),
  sig.level = 0.05,
  p.adj = FALSE,
  alternative = c("two.sided", "one.sided"),
  strict = TRUE
)
```

## Arguments

- object:

  a design object created in pwr4exp

- which:

  the factor of interest. Multiple factors can be combined using `:` or
  `*`, e.g., `"facA*facB"`, which represents a single factor that
  combines the levels of both factors.

- by:

  the variable to condition on

- contrast:

  A character string specifying the contrast method, one of "pairwise",
  "poly", or "trt.vs.ctrl". Alternatively, a numeric vector defining a
  single contrast or a (named) list of vectors specifying multiple
  custom contrasts. If a list is provided, each element must be a vector
  whose length matches the number of levels of the factor in each `by`
  group. In multi-factor scenarios, factor levels are combined and
  treated as a single factor.

- sig.level:

  significance level, default 0.05

- p.adj:

  whether the sig.level should be adjusted using the Bonferroni method,
  default FALSE

- alternative:

  one- or two-sided test. Can be abbreviated.

- strict:

  use strict interpretation in two-sided case

## Value

For each `by` condition, returns a data frame containing the contrast
value (`effect`), degrees of freedom (`df`), type I error rate
(`sig.level`), `power`, and the test direction (by `alternative`). When
multiple `by` conditions are present, the results are returned as a
list.

## Examples

``` r
rcbd <- designRCBD(
  treatments = c(2, 2),
  label = list(facA = c("1", "2"), facB = c("1", "2")),
  blocks = 12,
  formula = ~ facA*facB + (1|block),
  means = c(32, 35, 30, 37),
  vcomp = 4,
  sigma2 = 6
)

# If contrast is not specified, pairwise comparisons are conducted
pwr.contrast(rcbd, which = "facA") # Marginal effect of facA
#>               effect df sig.level     power alternative
#> facA1 - facA2     -5 33      0.05 0.9999995   two.sided
pwr.contrast(rcbd, which = "facA", by = "facB") # Conditional effect of facA within levels of facB
#> $`facB = 1`
#>               effect df sig.level     power alternative
#> facA1 - facA2     -3 33      0.05 0.8293757   two.sided
#> 
#> $`facB = 2`
#>               effect df sig.level     power alternative
#> facA1 - facA2     -7 33      0.05 0.9999993   two.sided
#> 

# Custom contrast vector, identical to pairwise comparison
pwr.contrast(rcbd, which = "facA", contrast = c(1, -1))
#>         effect df sig.level     power alternative
#> c(1,-1)     -5 33      0.05 0.9999995   two.sided
pwr.contrast(rcbd, which = "facA", by = "facB", contrast = c(1, -1))
#> $`facB = 1`
#>         effect df sig.level     power alternative
#> c(1,-1)     -3 33      0.05 0.8293757   two.sided
#> 
#> $`facB = 2`
#>         effect df sig.level     power alternative
#> c(1,-1)     -7 33      0.05 0.9999993   two.sided
#> 

# A single factor combining two factors
pwr.contrast(
  rcbd,
  which = "facA*facB",
  contrast = list(
    A1B1vs.A2B1 = c(1, -1, 0, 0),
    A1B1vs.A2B2 = c(1, 0, 0, -1)
  )
)
#>             effect df sig.level     power alternative
#> A1B1vs.A2B1     -3 33      0.05 0.8293757   two.sided
#> A1B1vs.A2B2     -5 33      0.05 0.9980739   two.sided
```
