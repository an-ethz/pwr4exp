# pwr4exp

`pwr4exp` provides a flexible framework for calculating statistical power across a range of experimental designs, with a particular focus on applications in animal science and related fields. The package supports approximate F-tests for general linear hypotheses in linear mixed models and t-tests for specific contrasts, utilizing the Satterthwaite method to approximate degrees of freedom.

The development version adds support for various correlation structures in repeated measures and spatial data. It also introduces options for creating input templates, making it easier to set up and perform power analyses. Note that sample size determination has been removed from this version.

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

Performing power analysis in `pwr4exp` involves two main steps:

1.  Create a design object.

2.  Pass the design object to power calculators.

### Example Usage

Here's an example of how you can generate a design and calculate power:

``` r
library(pwr4exp)

# Step 1. Create a design
## Define a completely randomized design with repeated measures
## by specifying it's data structure

n_subject = 6 # Subjects per treatment
n_trt = 3 # Number of treatments
n_hour = 8 # Number of repeated measures (time points)
trt = c("CON", "TRT1", "TRT2")

df.rep <- data.frame(
  subject = as.factor(rep(seq_len(n_trt*n_subject), each = n_hour)),
  hour = as.factor(rep(seq_len(n_hour), n_subject*n_trt)),
  trt = rep(trt, each = n_subject*n_hour)
)

## Generate a template for required input

mkdesign(formula = ~ trt*hour, data = df.rep)

## Create the design object

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

# Step 2: Calculate power
## Omnibus test
pwr.anova(design.rep)

## Contrast test (pairwise comparisons by hour)
pwr.contrast(design.rep, which = "trt", by = "hour", contrast = "pairwise")
```

## Learn More

To learn more about power analysis with `pwr4exp`, refer to the [vignette](https://an-ethz.github.io/pwr4exp/articles/pwr4exp.html) which contains:

-   Fundamental concepts of statistical power in linear mixed models.
-   Instructions for preparing and providing the required inputs.
-   Examples of power calculations for standard designs available in the package.
-   Examples of power analysis for more complex designs.

The package documentation is being updated. For any questions or suggestions, please feel free to contact the package maintainer.
