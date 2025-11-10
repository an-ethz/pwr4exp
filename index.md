# pwr4exp

`pwr4exp` supports statistical power calculations for diverse
experimental designs analyzed using linear mixed models.

It provides approximate **F-tests** for omnibus hypothesis tests and
**t-tests** for specific contrasts, employing the Satterthwaite method
for approximating degrees of freedom.

Various correlation structures (R-side) defined in the `nlme` package
can be used to model residuals.

## Installation

``` r
# You can install pwr4exp from CRAN
install.packages("pwr4exp")

# Or the development version from GitHub:
# install.packages("devtools")
devtools::install_github("an-ethz/pwr4exp")
```

``` r
# Load the library
library(pwr4exp)
```

## Functions

Performing a power analysis in `pwr4exp` involves two main steps:

### Step 1. **Define the design**

#### The `mkdesign` function

`mkdesign` is the most flexible function for defining an experimental
design based on its data structure, statistical model, treatment
effects, and variance-covariance components.

For example, the following code defines a completely randomized design
with repeated measures:

``` r
# Create a data frame reflecting the data structure
n_trt = 3 # Number of treatments
trt = c("CON", "TRT1", "TRT2")
n_subject = 6 # Number of subjects per treatment
n_hour = 8 # Number of repeated measures (time points)

df.rep <- data.frame(
  subject = as.factor(rep(seq_len(n_trt*n_subject), each = n_hour)),
  hour = as.factor(rep(seq_len(n_hour), n_subject*n_trt)),
  trt = rep(trt, each = n_subject*n_hour)
)

# Create the design object
crd.rep <- mkdesign(
  formula = ~ trt*hour, # model
  data = df.rep, # a data frame reflecting the data structure
  # treatment means at each hour reflecting effects
  means =  c(1, 2.50, 3.5,
             1, 3.50, 4.54,
             1, 3.98, 5.80,
             1, 4.03, 5.4,
             1, 3.68, 5.49,
             1, 3.35, 4.71,
             1, 3.02, 4.08,
             1, 2.94, 3.78), 
  sigma2 = 2, # residual variance
  # residual correlation structure, AR1
  correlation = corAR1(value = 0.6, form = ~ hour|subject) 
)
```

`pwr4exp` also provides functions for some common standard designs,
allowing a design to be defined by specifying key characteristics like
treatment structure and replications, without manually creating a data
frame.

``` r
# Create a completely randomized design
crd <- designCRD(
  treatments = 4, # one treatment factor with 4 levels
  replicates = 12, # 12 experimental units per treatment
  means = c(30, 28, 33, 35), # treatment means
  sigma2 = 10 # error variance
)

# Create a randomized complete block design
rcbd <- designRCBD(
  # two factors, each with 2 levels (2x2 factorial)
  treatments = c(2, 2), 
  blocks = 10, # 10 blocks
  means = c(30, 28, 33, 35), # treatment means (cell means)
  vcomp = 6, # block variance
  sigma2 = 4 # error variance
)

# Create a Latin Square design
lsd <- designLSD(
  # two factors, with 2 and 3 levels, respectively (2x3 factorial)
  treatments = c(2, 3), 
  squares = 4, # four squares
  reuse = "col", # column blocks are identical acorss squares
  means = c(30, 28, 33, 35, 34, 35), # treatment means (cell means)
  vcomp = c(5, 2), # variances of blocking factor (row and column)
  sigma2 = 3 # error variance
)

# Create a split-plot design
spd <- designSPD(
  trt.main = 2, # one factor with two levels at main plot
  trt.sub = 2, # one factor with two levels at subplot
  # 10 units per main plot factor (blocks at subplot level)
  replicates = 10, 
  means = c(30, 28, 33, 35), # treatment means (cell means)
  vcomp = 7, # main plot error
  sigma2 = 3 # residual error (subplot error)
)

# If not specified, these functions internally use a default model formula that includes main effects and all interactions, with block factors fitted as random effects.
```

### Step 2. **Calculate power**

Once the design has been correctly defined, the design object can be
passed to power calculation functions, including:

#### `pwr.anova`

Computes the power of F-tests for omnibus hypotheses.

``` r
pwr.anova(crd.rep)
```

#### `pwr.contrast`

Computes the power of t-tests for specific contrasts.

``` r
pwr.contrast(crd.rep, which = "trt", by = "hour", contrast = "trt.vs.ctrl")
```

## Learn More

To learn more about power analysis with `pwr4exp`, refer to the
[vignette](https://an-ethz.github.io/pwr4exp/articles/pwr4exp.html)
which contains:

- Instructions on preparing and providing the required inputs.
- Examples of power calculations for standard designs available in the
  package.
- Examples of power analysis for customized designs.
- Fundamental concepts of statistical power in linear mixed models.

For questions or suggestions, please open an
[issue](https://github.com/an-ethz/pwr4exp/issues) on our GitHub
repository or contact the package maintainer.
