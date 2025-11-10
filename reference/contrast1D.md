# Computes power of t-test for one-dimensional contrast matrices

Computes power of t-test for one-dimensional contrast matrices

## Usage

``` r
contrast1D(
  object,
  L,
  method = c("Satterthwaite"),
  sig.level = 0.05,
  alternative = c("two.sided", "one.sided"),
  strict = TRUE
)
```

## Arguments

- object:

  design object

- L:

  contrast vector

- method:

  DF approximation method, only "Satterthwaite" available currently

- sig.level:

  significance level, default 0.05

- alternative:

  one- or two-sided test

- strict:

  whether or not use strict interpretation in two-sided case

## Value

A data frame with columns for effect size, degrees of freedom,
significance level, power, and test type.
