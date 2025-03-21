% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/df_utils.R
\name{df.crd}
\alias{df.crd}
\title{Create a data frame of completely randomized design}
\usage{
df.crd(treatments, label, replicates)
}
\arguments{
\item{treatments}{An integer vector where each element represents the number of levels
of the corresponding treatment factor. A single integer (e.g., \code{treatments = n})
specifies one treatment factor with \code{n} levels. When multiple factors are provided,
they are arranged in a factorial treatment factor design. For example,
\code{treatments = c(2, 3)} creates a 2x3 factorial design with the first factor having 2 levels
and the second factor having 3 levels.}

\item{label}{Optional. A list of character vectors, each corresponding to a treatment factor.
The name of each vector specifies the factor's name, and its elements provide the labels for that factor's levels.
If no labels are provided, default labels will be used. For a single treatment factor, the default is
\code{list(trt = c("1", "2", ...))}, and for two treatment factors, the default is
\code{list(facA = c("1", "2", ...), facB = c("1", "2", ...))}.
For split-plot designs, the defaults are similar but include the ".main" and ".sub" suffixes for main plot and subplot factors.
For example:
\code{list(trt.main = c("1", "2", ...), trt.sub = c("1", "2", ...))}
\code{list(facA.main = c("1", "2", ...), facB.main = c("1", "2", ...),
      facA.sub = c("1", "2", ...), facB.sub = c("1", "2", ...))}
Label sets should be arranged so that the main plot factors come first, followed by the subplot factors.}

\item{replicates}{The number of experimental units per treatment.}
}
\value{
a data.frame representing the data structure of the design
}
\description{
Create a data frame of completely randomized design
}
