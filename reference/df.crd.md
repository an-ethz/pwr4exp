# Create a data frame of completely randomized design

Create a data frame of completely randomized design

## Usage

``` r
df.crd(treatments, label, replicates)
```

## Arguments

- treatments:

  An integer vector where each element represents the number of levels
  of the corresponding treatment factor. A single integer (e.g.,
  `treatments = n`) specifies one treatment factor with `n` levels. When
  multiple factors are provided, they are arranged in a factorial
  treatment factor design. For example, `treatments = c(2, 3)` creates a
  2x3 factorial design with the first factor having 2 levels and the
  second factor having 3 levels.

- label:

  Optional. A list of character vectors, each corresponding to a
  treatment factor. The name of each vector specifies the factor's name,
  and its elements provide the labels for that factor's levels. If no
  labels are provided, default labels will be used. For a single
  treatment factor, the default is `list(trt = c("1", "2", ...))`, and
  for two treatment factors, the default is
  `list(facA = c("1", "2", ...), facB = c("1", "2", ...))`. For
  split-plot designs, the defaults are similar but include the ".main"
  and ".sub" suffixes for main plot and subplot factors. For example:
  `list(trt.main = c("1", "2", ...), trt.sub = c("1", "2", ...))`
  `list(facA.main = c("1", "2", ...), facB.main = c("1", "2", ...), facA.sub = c("1", "2", ...), facB.sub = c("1", "2", ...))`
  Label sets should be arranged so that the main plot factors come
  first, followed by the subplot factors.

- replicates:

  The number of experimental units per treatment.

## Value

a data.frame representing the data structure of the design
