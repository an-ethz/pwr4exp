
# pwr4exp

`pwr4exp` is an R package that provides functions for power calculation and sample size determination in animal experiments. The package emphasizes the importance of specifying models for conducting power analyses and supports power analyses for main effects, interactions, and specific contrasts. Additionally, `pwr4exp` offers a flexible framework to perform power analysis on customized designs which are currently not predefined in the package.

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the development version of `pwr4exp` from [GitHub](https://github.com/WangKai7kkw/pwr4exp) with:

``` r
devtools::install_github("WangKai7kkw/pwr4exp", dependencies = TRUE, build_vignettes = TRUE)
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

## Learn More

To learn more about `pwr4exp`, read through the [vignette](https://github.com/WangKai7kkw/pwr4exp/vignettes/Introduction.Rmd) for `pwr4exp` which contains:

- Details on how to provide required inputs.
- Examples of calculating power and determining sample size for completely randomized design, random complete block design, Latin Square design.
- Examples using `pwr4exp` to assess the power of a customized design.

Additionally, more information is available at the [package website](wangkai7kkw.github.io/pwr4exp/).

The documentation for this package is being updated. If you have any questions or suggestions, please feel free to contact the package maintainer.



