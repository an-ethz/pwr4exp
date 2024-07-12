
# pwr4exp

`pwr4exp` is an R package that provides functions for power calculation and sample size determination in animal experiments. The package emphasizes the importance of specifying models for conducting power analyses and supports power analyses for main effects, interactions, and specific contrasts. Additionally, `pwr4exp` offers a flexible framework to perform power analysis on customized designs which are currently not predefined in the package.

<!-- badges: start -->
<!-- badges: end -->

pwr4exp is an R package that provides functions for power calculation and sample size determination  in animal experiments. pwr4exp emphasizes the importance of specifying models when conducting power analyses and enables power analyses for main effects, interactions, and specific contrasts. Additionally, pwr4exp offers a flexible framework that allows users to perform power analysis on customized designs which are not defined in the package curretly.

## Installation

You can install the development version of `pwr4exp` from [GitHub](https://github.com/yourgithub/pwr4exp) with:

``` r
devtools::install_github("WangKai7kkw/pwr4exp")
devtools::install_github(dependencies = TRUE, build_vignettes = TRUE)
```

## Functions

Performing power analysis in `pwr4exp` involves the following steps:
- First, create a desired design object (a list) using design generating functions.
- Once the design object is created, calculating power or determining sample size using `pwr4exp` is straightforward. Simply pass the design object to the power calculator for main effects and interactions, `pwr.anova()`, or for contrasts, `pwr.contrast()`.
- To determine the minimal sample size to achieve a target power, a quoted design object can be passed to the function `find_sample_size()`.

### Example Usage

Here's an example of how you can generate a design and calculate power:

```r
library(pwr4exp)
# Define a design
my_design <- designCRD(
  treatments = c(2, 2),
  replications = 10,
  beta = c(470, 30, -55, 5),
  sigma2 = 6400
)

# Calculate power for main effects and interaction
pwr.anova(design = my_design)

# Calculate power for contrasts
pwr.contrast(design = mydesign, specs = ~ facB | facA, contrast = "pairwise")
pwr.contrast(design = mydesign, specs = ~ facB * facA, contrast = "pairwise")
```

