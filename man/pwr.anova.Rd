% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/power-calculation.R
\name{pwr.anova}
\alias{pwr.anova}
\title{Power of omnibus tests}
\usage{
pwr.anova(object, sig.level = 0.05, type = c("III", "II", "I", "3", "2", "1"))
}
\arguments{
\item{object}{a design object created in pwr4exp}

\item{sig.level}{significance level, default 0.05}

\item{type}{the type of ANOVA table requested, default Type III}
}
\value{
a data frame with numerator degrees of freedom (NumDF), denominator
degrees of freedom (DenDF), type I error rate (sig.level), and power.
}
\description{
Calculate power for testing overall effects of treatment factors and their
interactions, i.e., statistical power of F-test.
}
\examples{
# generate an RCBD
rcbd <- designRCBD(
  treatments = c(2, 2),
  label = list(facA = c("1", "2"), facB = c("1", "2")),
  blocks = 12,
  formula = ~ facA*facB + (1|block),
  means = c(32, 35, 30, 37),
  vcomp = 4,
  sigma2 = 6
)
# power of omnibus test
pwr.anova(rcbd)

}
\seealso{
\code{\link[=mkdesign]{mkdesign()}}, \code{\link[=designCRD]{designCRD()}}, \code{\link[=designRCBD]{designRCBD()}}, \code{\link[=designLSD]{designLSD()}}, \code{\link[=designCOD]{designCOD()}}, \code{\link[=designSPD]{designSPD()}}, \code{\link[=pwr.summary]{pwr.summary()}} and \code{\link[=pwr.contrast]{pwr.contrast()}}
}
