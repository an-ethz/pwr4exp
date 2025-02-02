% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power-calculation.R
\name{pwr.contrast}
\alias{pwr.contrast}
\title{Power of contrasts}
\usage{
pwr.contrast(
  object,
  which,
  by = NULL,
  contrast = c("pairwise", "poly", "trt.vs.ctrl"),
  sig.level = 0.05,
  p.adj = FALSE,
  alternative = c("two.sided", "one.sided"),
  strict = TRUE
)
}
\arguments{
\item{object}{design object created in pwr4exp}

\item{which}{the factor of interest}

\item{by}{the variable to condition on}

\item{contrast}{contrast method, include "pairwise", "poly", and "trt.vs.ctrl",
or any manually defined contrast vector}

\item{sig.level}{significance level, default 0.05}

\item{p.adj}{whether the sig.level should be adjusted using the Bonferroni method, default FALSE}

\item{alternative}{one- or two-sided test. Can be abbreviated.}

\item{strict}{use strict interpretation in two-sided case}
}
\value{
a data frame or a list of data frame separated by conditions.
}
\description{
Calculate power for testing various contrasts.
}
\examples{
rcbd <- designRCBD(
  treatments = c(2, 2),
  label = list(facA = c("1", "2"), facB = c("1", "2")),
  blocks = 12,
  formula = ~ facA*facB + (1|block),
  means = c(32, 35, 30, 37),
  vcomp = 4,
  sigma2 = 6
)
pwr.contrast(rcbd, which = "facA", by = "facB")
}
