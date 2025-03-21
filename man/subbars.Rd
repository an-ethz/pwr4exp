% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{subbars}
\alias{subbars}
\title{Substitute Bars}
\usage{
subbars(term)
}
\arguments{
\item{term}{a mixed-model formula}
}
\value{
the formula with all |  and || operators replaced by +
}
\description{
Substitute the '+' function for the '|' and '||' function in a mixed-model
formula.
}
\section{Note}{
 This function is called recursively on individual
terms in the model, which is why the argument is called \code{term} and not
a name like \code{form}, indicating a formula.
}

\keyword{functions}
\keyword{internal}
\keyword{lme4}
