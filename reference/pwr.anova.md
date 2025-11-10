# Power of omnibus tests

Calculates the statistical power for testing the overall effects of
treatment factors and their interactions, i.e., power of F-test.

## Usage

``` r
pwr.anova(object, sig.level = 0.05, type = c("III", "II", "I", "3", "2", "1"))
```

## Arguments

- object:

  a design object created in pwr4exp

- sig.level:

  significance level, default 0.05

- type:

  the type of ANOVA table requested, default Type III

## Value

a data frame with numerator degrees of freedom (`NumDF`), denominator
degrees of freedom (`DenDF`), type I error rate (`sig.level`), and
`power`.

## See also

[mkdesign](https://an-ethz.github.io/pwr4exp/reference/mkdesign.md),
[designCRD](https://an-ethz.github.io/pwr4exp/reference/create_designs.md),
[designRCBD](https://an-ethz.github.io/pwr4exp/reference/create_designs.md),
[designLSD](https://an-ethz.github.io/pwr4exp/reference/create_designs.md),
[designCOD](https://an-ethz.github.io/pwr4exp/reference/create_designs.md),
[designSPD](https://an-ethz.github.io/pwr4exp/reference/create_designs.md),
[pwr.summary](https://an-ethz.github.io/pwr4exp/reference/pwr.summary.md)
and
[pwr.contrast](https://an-ethz.github.io/pwr4exp/reference/pwr.contrast.md)

## Examples

``` r
# generate an RCBD
rcbd <- designRCBD(
  treatments = c(2, 2),
  label = list(facA = c("1", "2"), facB = c("1", "2")),
  blocks = 12,
  formula = ~ facA*facB + (1|block),
  means = c(32, 35, 30, 37),
  vcomp = 4,
  sigma2 = 6
)
# power of omnibus test
pwr.anova(rcbd)
#> Power of type III F-test
#>           NumDF DenDF sig.level   power
#> facA          1    33      0.05 1.00000
#> facB          1    33      0.05 0.05000
#> facA:facB     1    33      0.05 0.78387
```
