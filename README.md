# pwr4exp

The `pwr4exp` R package provides functions for power calculation and sample size determination in standard experimental designs in animal science and beyond. The package emphasizes the importance of specifying models for conducting power analyses and supports power analyses for main effects, interactions, and specific contrasts. Additionally, `pwr4exp` offers a flexible framework to perform power analysis on customized designs which are currently not predefined in the package. The development version incorporates various correlation structures for repeated measures and spatial data directly from the `nlme` package.

<!-- badges: start -->

<!-- badges: end -->

## Installation

``` r
# You can install pwr4exp from CRAN
install.packages("pwr4exp")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("an-ethz/pwr4exp")
```

## Functions

Performing power analysis in `pwr4exp` involves the following steps: - First, create a design object using the design generating functions. - Once the design object is created, calculating power or determining sample size using `pwr4exp` is straightforward. Simply pass the design object to the power calculator for main effects and interactions, `pwr.anova()`, or for contrasts, `pwr.contrast()`.

### Example Usage

Here's an example of how you can generate a design and calculate power:

``` r
library(pwr4exp)

## representing a CRD with repeated measures at 8 time points

n_subject = 6
n_trt = 3
n_hour = 8
trt = c("CON", "TRT1", "TRT2")

df.rep <- data.frame(
  subject = as.factor(rep(seq_len(n_trt*n_subject), each = n_hour)),
  hour = as.factor(rep(seq_len(n_hour), n_subject*n_trt)),
  trt = rep(trt, each = n_subject*n_hour)
)

# templates
mkdesign(formula = ~ trt*hour, data = df.rep)

# create the design

design.rep <- mkdesign(
formula = ~ trt*hour,
data = df.rep,
means =  c(1, 2.50, 3.5,
           1, 3.50, 4.54,
           1, 3.98, 5.80,
           1, 4.03, 5.4,
           1, 3.68, 5.49,
           1, 3.35, 4.71,
           1, 3.02, 4.08,
           1, 2.94, 3.78),
sigma2 = 2,
correlation = corAR1(value = 0.6, form = ~ hour|subject)
)

pwr.anova(design.rep)
```

## Learn More

The page will be updated recently to reflect the newly add functions.

~~To learn more about `pwr4exp`, read through the [vignette](https://an-ethz.github.io/pwr4exp/articles/pwr4exp.html) for `pwr4exp` which contains:~~

~~- Foundational Concepts.~~ ~~- Instructions for providing the required inputs.~~ \~\~- Examples of conducting power analysis for the standard designs available in the package. ~~- Examples of using `pwr4exp` to assess the power of customized designs.~~

The documentation for this package is being updated. If you have any questions or suggestions, please feel free to contact the package maintainer.
