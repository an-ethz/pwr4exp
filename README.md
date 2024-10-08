
# pwr4exp

The `pwr4exp` R package provides functions for power calculation and sample size determination in standard experimental designs in animal science and beyond. The package emphasizes the importance of specifying models for conducting power analyses and supports power analyses for main effects, interactions, and specific contrasts. Additionally, `pwr4exp` offers a flexible framework to perform power analysis on customized designs which are currently not predefined in the package.

<!-- badges: start -->
<!-- badges: end -->

## Installation

```r
# You can install pwr4exp from CRAN
install.packages("pwr4exp")

# Or the the development version from GitHub:
# install.packages("devtools")
devtools::install_github("an-ethz/pwr4exp")
```

## Functions

Performing power analysis in `pwr4exp` involves the following steps:
- First, create a design object using the design generating functions.
- Once the design object is created, calculating power or determining sample size using `pwr4exp` is straightforward. Simply pass the design object to the power calculator for main effects and interactions, `pwr.anova()`, or for contrasts, `pwr.contrast()`.
- To determine the minimal sample size to achieve a target power, a quoted design object can be passed to the function `find_sample_size()`.

### Example Usage

Here's an example of how you can generate a design and calculate power:

```r
library(pwr4exp)
# Create a design
rcbd <- designRCBD(
  treatments = c(2, 2),
  blocks = 10,
  beta = c(470, 30, -55, 5),
  VarCov = 3200,
  sigma2 = 3200
)
# Calculate power for the ominubus test (i.e., F-test)
pwr.anova(design = rcbd)

# Calculate power for specific contrasts
pwr.contrast(design = rcbd, specs = ~ facB | facA, contrast = "pairwise")
```

## Learn More

To learn more about `pwr4exp`, read through the [vignette](https://an-ethz.github.io/pwr4exp/articles/pwr4exp.html) for `pwr4exp` which contains:

- Foundational Concepts.
- Instructions for providing the required inputs.
- Examples of conducting power analysis for the standard designs available in the package.
- Examples of using `pwr4exp` to assess the power of customized designs.

The documentation for this package is being updated. If you have any questions or suggestions, please feel free to contact the package maintainer.
