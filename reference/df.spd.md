# Create data frame for split-plot design

Create data frame for split-plot design

## Usage

``` r
df.spd(trt.main, trt.sub, label, replicates)
```

## Arguments

- trt.main:

  an integer-valued vector specifying the treatment structure at main
  plot level, similar to
  [`df.crd`](https://an-ethz.github.io/pwr4exp/reference/df.crd.md).

- trt.sub:

  an integer-valued vector specifying the treatment structure at sub
  plot level, similar to `trt.main`.

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
  `list(trt.main = c("1", "2", ...), trt.sub = c("1", "2", ...))` and
  `list(facA.main = c("1", "2", ...), facB.main = c("1", "2", ...), facA.sub = c("1", "2", ...), facB.sub = c("1", "2", ...))`.
  Label sets should be arranged so that the main plot factors come
  first, followed by the subplot factors.

- replicates:

  the number of experimental units (main plots) per treatment of main
  plot factors.

## Value

a data.frame representing the data structure of the design
