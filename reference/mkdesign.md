# Create a Design Object for Power Calculation

Generate a design object for power analysis by specifying a model
formula and data frame. This object is not a true experimental design as
created by design generation procedures, where randomization and unit
allocation are performed. Instead, it serves as an object containing all
necessary information for power analysis, including design matrices,
assumed values of model effects, and other internally calculated
parameters.

## Usage

``` r
mkdesign(
  formula,
  data,
  beta = NULL,
  means = NULL,
  vcomp = NULL,
  sigma2 = NULL,
  correlation = NULL,
  template = FALSE,
  REML = TRUE
)
```

## Arguments

- formula:

  A right-hand-side [formula](https://rdrr.io/r/stats/formula.html)
  specifying the model for testing treatment effects, with terms on the
  right of [~](https://rdrr.io/r/base/tilde.html) , following lme4::lmer
  syntax for random effects.

- data:

  A data frame with all independent variables specified in the model,
  matching the design's structure.

- beta:

  One of the optional inputs for fixed effects. A vector of model
  coefficients where factor variable coefficients correspond to dummy
  variables created using treatment contrast (stats::contr.treatment).

- means:

  One of the optional inputs for fixed effects. A vector of marginal or
  conditioned means (if factors have interactions). Regression
  coefficients are required for numerical variables. Either `beta` or
  `means` must be provided, and their values must strictly follow a
  specific order. A template can be created to indicate the required
  input values and their order. See "Details" for more information.

- vcomp:

  A vector of variance-covariance components for random effects, if
  present. The values must follow a strict order. See "Details".

- sigma2:

  error variance.

- correlation:

  Specifies residual (R-side) correlation structures using
  [nlme::corClasses](https://rdrr.io/pkg/nlme/man/corClasses.html)
  functions. See "Details" for more information.

- template:

  Default is `FALSE`. If `TRUE` or when only the formula and data are
  provided, a template for `beta`, `means`, and `vcomp` is generated to
  indicate the required input order.

- REML:

  Specifies whether to use REML or ML information matrix. Default is
  `TRUE`.

## Value

A list object containing all essential components for power calculation.
This includes:

- Structural components (deStruct): including design matrices for fixed
  and random effects, variance-covariance matrices for random effects
  and residuals, etc.

- Internally calculated higher-level parameters (deParam), including
  variance-covariance matrix of beta coefficients (vcov_beta),
  variance-covariance matrix of variance parameters (vcov_varpar),
  gradient matrices (Jac_list), etc.

## Details

- **data**: A long-format data frame is required, as typically used in R
  for fitting linear models. This data frame can be created manually or
  with the help of design creation packages such as agricolae,
  AlgDesign, crossdes, or FrF2. It should include all independent
  variables specified in the model (e.g., treatments, blocks, subjects).
  All the irrelevant variables not specified in the model are ignored.

- **template**: Templates are automatically generated when only the
  formula and data are supplied, or explicitly if `template = TRUE`.
  Templates serve as guides for specifying inputs:

  - **Template for `beta`**: Represents the sequence of model
    coefficients.

  - **Template for `means`**: Specifies the order of means (for
    categorical variables) and/or regression coefficients (for
    continuous variables), depending on the scenario:

    - *Categorical variables without interactions*: Requires marginal
      means for each level of the categorical variable(s).

    - *Interactions among categorical variables*: Requires conditional
      (cell) means for all level combinations.

    - *Numerical variables without interactions*: Requires regression
      coefficients. The intercept must also be included if there are no
      categorical variables in the model.

    - *Interactions among numerical variables*: Requires regression
      coefficients for both main effects and interaction terms. The
      intercept must also be included if there are no categorical
      variables in the model.

    - *Categorical-by-numerical interactions*: Requires regression
      coefficients for the numerical variable at each level of the
      categorical variable, as well as marginal means for the levels of
      the categorical variable.

    Note: For models containing only numerical variables, the inputs for
    `means` and `beta` are identical. See the "Examples" for
    illustrative scenarios.

  - **Template for `vcomp`**: Represents a variance-covariance matrix,
    where integers indicate the order of variance components in the
    input vector.

- **correlation**: Various residual correlation structures can be
  specified following the instructions from
  [corClasses](https://an-ethz.github.io/pwr4exp/reference/corClasses.md)
  constructors.

  Note: In the original
  [`corAR1`](https://rdrr.io/pkg/nlme/man/corAR1.html),
  [`corARMA`](https://rdrr.io/pkg/nlme/man/corARMA.html), and
  [`corSymm`](https://rdrr.io/pkg/nlme/man/corSymm.html) functions, the
  covariate `t` in the correlation formula `~ t` or `~ t | g` must be an
  integer class. In pwr4exp, the covariate can also be a factor class,
  which is then converted to an integer internally for sorting purposes.
  The class of the covariate variable in the model formula, if present,
  will not be converted. For example, a time covariate can be fitted as
  a factor in the model formula, whereas it is converted to an integer
  in the correlation formula temporarily for matrix sorting.

## Examples

``` r
# Using templates for specifying "means"

# Create an example data frame with four categorical variables (factors)
# and two numerical variables
df1 <- expand.grid(
  fA = factor(1:2),
  fB = factor(1:2),
  fC = factor(1:3),
  fD = factor(1:3),
  subject = factor(1:10)
)
df1$x <- rnorm(nrow(df1))  # Numerical variable x
df1$z <- rnorm(nrow(df1))  # Numerical variable z

## Categorical variables without interactions
# Means of each level of fA and fB are required in sequence.
mkdesign(~ fA + fB, df1)$fixeff$means
#> fA1 fA2 fB1 fB2 
#>   1   2   3   4 

## Interactions among categorical variables
# Cell means for all combinations of levels of fA and fB are required.
mkdesign(~ fA * fB, df1)$fixeff$means
#> fA1:fB1 fA2:fB1 fA1:fB2 fA2:fB2 
#>       1       2       3       4 

## Numerical variables without and with interactions, identical to beta.
# Without interactions: Regression coefficients are required
mkdesign(~ x + z, df1)$fixeff$means
#> (Intercept)           x           z 
#>           1           2           3 

# With interactions: Coefficients for main effects and interaction terms are required.
mkdesign(~ x * z, df1)$fixeff$means
#> (Intercept)           x           z         x:z 
#>           1           2           3           4 

## Categorical-by-numerical interactions
# Marginal means for each level of fA, and regression coefficients for x
# at each level of fA are required.
mkdesign(~ fA * x, df1)$fixeff$means
#>   fA1   fA2 fA1:x fA2:x 
#>     1     2     3     4 

## Three factors with interactions:
# Cell means for all 12 combinations of the levels of fA, fB, and fC are required.
mkdesign(~ fA * fB * fC, df1)
#> $fixeff
#> $fixeff$beta
#> (Intercept)         fA2         fB2         fC2         fC3     fA2:fB2 
#>           1           2           3           4           5           6 
#>     fA2:fC2     fA2:fC3     fB2:fC2     fB2:fC3 fA2:fB2:fC2 fA2:fB2:fC3 
#>           7           8           9          10          11          12 
#> 
#> $fixeff$means
#> fA1:fB1:fC1 fA2:fB1:fC1 fA1:fB2:fC1 fA2:fB2:fC1 fA1:fB1:fC2 fA2:fB1:fC2 
#>           1           2           3           4           5           6 
#> fA1:fB2:fC2 fA2:fB2:fC2 fA1:fB1:fC3 fA2:fB1:fC3 fA1:fB2:fC3 fA2:fB2:fC3 
#>           7           8           9          10          11          12 
#> 
#> 
#> $varcov
#> NULL
#> 

# A design that mixes the above-mentioned scenarios:
# - Interactions among three categorical variables (fA, fB, fC)
# - A categorical-by-numerical interaction (fD * x)
# - Main effects for another numerical variable (z)
# The required inputs are:
# - Cell means for all combinations of levels of fA, fB, and fC
# - Means for each level of fD
# - Regression coefficients for x at each level of fD
# - Regression coefficients for z
mkdesign(~ fA * fB * fC + fD * x + z, df1)$fixeff$means
#>         fD1         fD2         fD3           z       fD1:x       fD2:x 
#>           1           2           3           4           5           6 
#>       fD3:x fA1:fB1:fC1 fA2:fB1:fC1 fA1:fB2:fC1 fA2:fB2:fC1 fA1:fB1:fC2 
#>           7           8           9          10          11          12 
#> fA2:fB1:fC2 fA1:fB2:fC2 fA2:fB2:fC2 fA1:fB1:fC3 fA2:fB1:fC3 fA1:fB2:fC3 
#>          13          14          15          16          17          18 
#> fA2:fB2:fC3 
#>          19 

# Using templates for specifying "vcomp"

# Assume df1 represents an RCBD with "subject" as a random blocking factor.
## Variance of the random effect "subject" (intercept) is required.
mkdesign(~ fA * fB * fC * fD + (1 | subject), df1)$varcov
#> $subject
#>             (Intercept)
#> (Intercept)           1
#> 

# Demonstration of templates for more complex random effects
## Note: This example is a demo and statistically incorrect for this data
## (no replicates under subject*fA). It only illustrates variance-covariance templates.
## Inputs required:
## - Variance of the random intercept (1st)
## - Covariance between the intercept and "fA2" (2nd)
## - Variance of "fA2" (3rd)
mkdesign(~ fA * fB * fC * fD + (1 + fA | subject), df1)$varcov
#> $subject
#>             (Intercept) fA2
#> (Intercept)           1   2
#> fA2                   2   3
#> 

# Power analysis for repeated measures

## Create a data frame for a CRD with repeated measures
n_subject <- 6
n_trt <- 3
n_hour <- 8
trt <- c("CON", "TRT1", "TRT2")
df2 <- data.frame(
  subject = as.factor(rep(seq_len(n_trt * n_subject), each = n_hour)), # Subject as a factor
  hour = as.factor(rep(seq_len(n_hour), n_subject * n_trt)),           # Hour as a factor
  trt = rep(trt, each = n_subject * n_hour)                           # Treatment as a factor
)

## Templates
temp <- mkdesign(formula = ~ trt * hour, data = df2)
temp$fixeff$means  # Fixed effects means template
#>  trtCON:hour1 trtTRT1:hour1 trtTRT2:hour1  trtCON:hour2 trtTRT1:hour2 
#>             1             2             3             4             5 
#> trtTRT2:hour2  trtCON:hour3 trtTRT1:hour3 trtTRT2:hour3  trtCON:hour4 
#>             6             7             8             9            10 
#> trtTRT1:hour4 trtTRT2:hour4  trtCON:hour5 trtTRT1:hour5 trtTRT2:hour5 
#>            11            12            13            14            15 
#>  trtCON:hour6 trtTRT1:hour6 trtTRT2:hour6  trtCON:hour7 trtTRT1:hour7 
#>            16            17            18            19            20 
#> trtTRT2:hour7  trtCON:hour8 trtTRT1:hour8 trtTRT2:hour8 
#>            21            22            23            24 

## Create a design object
# Assume repeated measures within a subject follow an AR1 process with a correlation of 0.6
design <- mkdesign(
  formula = ~ trt * hour,
  data = df2,
  means = c(1, 2.50, 3.50,
            1, 3.50, 4.54,
            1, 3.98, 5.80,
            1, 4.03, 5.40,
            1, 3.68, 5.49,
            1, 3.35, 4.71,
            1, 3.02, 4.08,
            1, 2.94, 3.78),
  sigma2 = 2,
  correlation = corAR1(value = 0.6, form = ~ hour | subject)
)

pwr.anova(design)  # Perform power analysis
#> Power of type III F-test
#>          NumDF  DenDF sig.level   power
#> trt          2 21.563      0.05 1.00000
#> hour         7 86.055      0.05 0.74687
#> trt:hour    14 86.055      0.05 0.38500

## When time is treated as a numeric variable
# Means of treatments and regression coefficients for hour at each treatment level are required
df2$hour <- as.integer(df2$hour)
mkdesign(formula = ~ trt * hour, data = df2)$fixeff$means
#>       trtCON      trtTRT1      trtTRT2  trtCON:hour trtTRT1:hour trtTRT2:hour 
#>            1            2            3            4            5            6 

## Polynomial terms of time in the model
mkdesign(formula = ~ trt + hour + I(hour^2) + trt:hour + trt:I(hour^2), data = df2)$fixeff$means
#>            trtCON           trtTRT1           trtTRT2       trtCON:hour 
#>                 1                 2                 3                 4 
#>      trtTRT1:hour      trtTRT2:hour  trtCON:I(hour^2) trtTRT1:I(hour^2) 
#>                 5                 6                 7                 8 
#> trtTRT2:I(hour^2) 
#>                 9 
```
