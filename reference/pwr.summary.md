# Power for model coefficients

Computes the statistical power for testing (t-test) model coefficients.

## Usage

``` r
pwr.summary(object, sig.level = 0.05)
```

## Arguments

- object:

  design object

- sig.level:

  significance level, default 0.05

## Value

a data frame containing model coefficients, degrees of freedom (`df`),
type I error rate (`sig.level`), power, and the test direction
(`alternative`).

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
pwr.summary(rcbd)
#>             effect       df sig.level     power alternative
#> (Intercept)     32 29.72973      0.05 1.0000000   two.sided
#> facA2            3 33.00000      0.05 0.8293757   two.sided
#> facB2           -2 33.00000      0.05 0.4927485   two.sided
#> facA2:facB2      4 33.00000      0.05 0.7838664   two.sided
```
