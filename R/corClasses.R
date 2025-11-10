###############################################################################
# This file is part of the pwr4exp package.
#
# pwr4exp is licensed under the GNU General Public License (GPL);
# see the LICENSE file for details.
#
# pwr4exp is provided in the hope that it will be useful, but WITHOUT ANY WARRANTY,
# including the implied warranties of MERCHANTABILITY and FITNESS FOR A PARTICULAR PURPOSE.
# For more information, please refer to the GNU General Public License.
###############################################################################

#' Correlation Structure Classes
#'
#' @description
#' Standard classes of correlation structures (\code{corStruct}) available from the \pkg{nlme} package
#' and re-exported by \pkg{pwr4exp} for convenience when specifying correlation structures in \code{\link{mkdesign}}.
#'
#' All arguments are identical to the corresponding \pkg{nlme} functions. For more details on the original implementations,
#' see [nlme::corClasses].
#'
#' Note: In the original [nlme::corAR1], [nlme::corARMA], and [nlme::corSymm] functions,
#' the covariate \code{t} in the correlation formula \code{~ t} or \code{~ t | g} must be an integer class. In \pkg{pwr4exp},
#' the covariate can also be a factor class, which is then converted to an integer internally for sorting purposes.
#' The class of the covariate variable in the model formula, if present, will not be converted. For example, a time covariate
#' can be fitted as a factor in the model formula, whereas it is converted to an integer in the correlation formula temporarily for matrix sorting.
#'
#' @details
#' Available standard classes:
#' \describe{
#'   \item{\code{\link{corAR1}}}{autoregressive process of order 1.}
#'   \item{\code{\link{corARMA}}}{autoregressive moving average process, with arbitrary orders
#'     for the autoregressive and moving average components.}
#'   \item{\code{\link{corCAR1}}}{continuous AR(1)}
#'   \item{\code{\link{corCompSymm}}}{compound symmetry structure corresponding to a constant correlation.}
#'   \item{\code{\link{corExp}}}{exponential spatial correlation.}
#'   \item{\code{\link{corGaus}}}{Gaussian spatial correlation.}
#'   \item{\code{\link{corLin}}}{linear spatial correlation.}
#'   \item{\code{\link{corRatio}}}{Rational quadratics spatial correlation.}
#'   \item{\code{\link{corSpher}}}{spherical spatial correlation.}
#'   \item{\code{\link{corSymm}}}{general correlation matrix, with no additional structure.}
#' }
#'
#' @references
#' Pinheiro, J.C., and Bates, D.M. (2000) "Mixed-Effects Models in S and S-PLUS", Springer.
#'
#' @name corClasses
#' @rdname corClasses
NULL

#' @name corAR1
#' @rdname corAR1
#' @title AR(1) Correlation Structure
#' @description Re-exports [nlme::corAR1] from the \pkg{nlme} package.
#' @usage corAR1(value, form, fixed)
#' @param value numeric value for the correlation parameter.
#' @param form a one-sided formula of the form \code{~ t}, \code{~ 1 | g}, or \code{~ t | g},
#' specifying a time covariate \code{t}, a grouping factor \code{g}, or both.
#' When no time covariate is specified, the row order of the data within each group is assumed to
#' represent the chronological order of measurements.
#' @param fixed unused
#' @details
#' In the original [nlme::corAR1] function, a covariate \code{t} must be an integer class.
#' In \pkg{pwr4exp}, \code{t} can also be a factor class, which is then converted to an integer internally for chronological order.
#' The class of \code{t} in the model formula, if present, is not converted. For example, a time covariate can be fitted as a factor
#' in the model formula, whereas it is converted to an integer in the correlation formula.
#'
#' @export
#' @importFrom nlme corAR1
#' @seealso [corClasses]
#' See [nlme::corAR1] for original documentation.
corAR1 <- nlme::corAR1

#' @name corARMA
#' @rdname corARMA
#' @title ARMA(p,q) Correlation Structure
#' @description Re-exports [nlme::corARMA] from the \pkg{nlme} package.
#' @usage corARMA(value, form, p, q, fixed)
#' @param value numeric vector of parameter values
#' @param form a one-sided formula of the form \code{~ t}, \code{~ 1 | g}, or \code{~ t | g},
#' specifying a time covariate \code{t}, a grouping factor \code{g}, or both.
#' When no time covariate is specified, the row order of the data within each group is assumed to
#' represent the chronological order of measurements.
#' @param p non-negative integer specifying the autoregressive order
#' @param q non-negative integer specifying the moving average order
#' @param fixed unused
#' @details
#' In the original [nlme::corARMA] function, a covariate \code{t} must be an integer class.
#' In \pkg{pwr4exp}, \code{t} can also be a factor class, which is then converted to an integer internally for chronological order.
#' The class of \code{t} in the model formula, if present, is not converted. For example, a time covariate can be fitted as a factor
#' in the model formula, whereas it is converted to an integer in the correlation formula.
#' @export
#' @importFrom nlme corARMA
#' @seealso [corClasses]
#' See [nlme::corARMA] for original documentation.
corARMA <- nlme::corARMA

#' @name corCAR1
#' @rdname corCAR1
#' @title Continuous AR(1) Correlation Structure
#' @description Re-exports [nlme::corCAR1] from the \pkg{nlme} package.
#' @usage corCAR1(value, form, fixed)
#' @param value numeric value for the correlation parameter
#' @param form a one-sided formula of the form \code{~ t}, \code{~ 1 | g}, or \code{~ t | g},
#' specifying a time covariate \code{t}, a grouping factor \code{g}, or both.
#' @param fixed unused
#' @export
#' @importFrom nlme corCAR1
#' @seealso [corClasses].
#' See [nlme::corCAR1] for original documentation.
corCAR1 <- nlme::corCAR1

#' @name corCompSymm
#' @rdname corCompSymm
#' @title Compound Symmetry Correlation Structure
#' @description Re-exports [nlme::corCompSymm] from the \pkg{nlme} package.
#' @usage corCompSymm(value, form, fixed)
#' @param value numeric vector of parameter values
#' @param form a one-sided formula of the form \code{~ 1 | g} specifying a grouping factor.
#' @param fixed unused
#' @export
#' @importFrom nlme corCompSymm
#' @seealso [corClasses]
#' See [nlme::corCompSymm] for original documentation.
corCompSymm <- nlme::corCompSymm

#' @name corExp
#' @rdname corExp
#' @title Exponential Correlation Structure
#' @description Re-exports [nlme::corExp] from the \pkg{nlme} package.
#' @usage corExp(value, form, nugget, metric, fixed)
#' @param value numeric value(s) for the parameter(s)
#' @param form a one-sided formula of the form \code{~ S1 + ... + Sp}, or \code{~ S1 + ... + Sp | g}
#' spatial covariates S1 through Sp and, optionally, a grouping factor g.
#' @param nugget logical; if \code{TRUE} a nugget effect is added
#' @param metric character string specifying the distance metric
#' @param fixed unused
#' @export
#' @importFrom nlme corExp
#' @seealso [corClasses]
#' See [nlme::corExp] for original documentation.
corExp <- nlme::corExp

#' @name corGaus
#' @rdname corGaus
#' @title Gaussian Correlation Structure
#' @description Re-exports [nlme::corGaus] from the \pkg{nlme} package.
#' @usage corGaus(value, form, nugget, metric, fixed)
#' @param value numeric value(s) for the parameter(s)
#' @param form a one-sided formula of the form \code{~ S1 + ... + Sp}, or \code{~ S1 + ... + Sp | g}
#' spatial covariates S1 through Sp and, optionally, a grouping factor g.
#' @param nugget logical; if \code{TRUE} a nugget effect is added
#' @param metric character string specifying the distance metric
#' @param fixed unused
#' @export
#' @importFrom nlme corGaus
#' @seealso [corClasses]
#' See [nlme::corGaus] for original documentation.
corGaus <- nlme::corGaus

#' @name corLin
#' @rdname corLin
#' @title Linear Correlation Structure
#' @description Re-exports [nlme::corLin] from the \pkg{nlme} package.
#' @usage corLin(value, form, nugget, metric, fixed)
#' @param value numeric value(s) for the parameter(s)
#' @param form a one-sided formula of the form \code{~ S1 + ... + Sp}, or \code{~ S1 + ... + Sp | g}
#' spatial covariates S1 through Sp and, optionally, a grouping factor g.
#' @param nugget logical; if \code{TRUE} a nugget effect is added
#' @param metric character string specifying the distance metric
#' @param fixed unused
#' @export
#' @importFrom nlme corLin
#' @seealso [corClasses]
#' See [nlme::corLin] for original documentation.
corLin <- nlme::corLin

#' @name corRatio
#' @rdname corRatio
#' @title Rational Quadratics Correlation Structure
#' @description Re-exports [nlme::corRatio] from the \pkg{nlme} package.
#' @usage corRatio(value, form, nugget, metric, fixed)
#' @param value numeric value(s) for the parameter(s)
#' @param form a one-sided formula of the form \code{~ S1 + ... + Sp}, or \code{~ S1 + ... + Sp | g}
#' spatial covariates S1 through Sp and, optionally, a grouping factor g.
#' @param nugget logical; if \code{TRUE} a nugget effect is added
#' @param metric character string specifying the distance metric
#' @param fixed unused
#' @export
#' @importFrom nlme corRatio
#' @seealso [corClasses]
#' See [nlme::corRatio] for original documentation.
corRatio <- nlme::corRatio

#' @name corSpher
#' @rdname corSpher
#' @title Spherical Correlation Structure
#' @description Re-exports [nlme::corSpher] from the \pkg{nlme} package.
#' @usage corSpher(value, form, nugget, metric, fixed)
#' @param value numeric value(s) for the parameter(s)
#' @param form a one-sided formula of the form \code{~ S1 + ... + Sp}, or \code{~ S1 + ... + Sp | g}
#' spatial covariates S1 through Sp and, optionally, a grouping factor g.
#' @param nugget logical; if \code{TRUE} a nugget effect is added
#' @param metric character string specifying the distance metric
#' @param fixed unused
#' @export
#' @importFrom nlme corSpher
#' @seealso [corClasses]
#' See [nlme::corSpher] for original documentation.
corSpher <- nlme::corSpher

#' @name corSymm
#' @rdname corSymm
#' @title General Correlation Structure
#' @description Re-exports [nlme::corSymm] from the \pkg{nlme} package.
#' @usage corSymm(value, form, fixed)
#' @param value numeric vector of correlation parameter values
#' @param form a one-sided formula of the form \code{~ t}, \code{~ 1 | g}, or \code{~ t | g},
#' specifying a covariate \code{t}, a grouping factor \code{g}, or both. A covariate indicates
#' the order of the data rows in the residual error matrix of each group.
#' @param fixed unused
#' @export
#' @importFrom nlme corSymm
#' @seealso [corClasses]
#' See [nlme::corSymm] for original documentation.
corSymm <- nlme::corSymm
