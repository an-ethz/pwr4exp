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
#' These functions re-export \code{\link[nlme]{corClasses}} constructors from the \pkg{nlme} package.
#'
#' @name corClasses
#' @rdname corClasses
#' @export
#' @importFrom nlme corAR1
corAR1 <- nlme::corAR1

#' @export
#' @rdname corClasses
#' @importFrom nlme corARMA
corARMA <- nlme::corARMA

#' @export
#' @rdname corClasses
#' @importFrom nlme corCAR1
corCAR1 <- nlme::corCAR1

#' @export
#' @rdname corClasses
#' @importFrom nlme corCompSymm
corCompSymm <- nlme::corCompSymm

#' @export
#' @rdname corClasses
#' @importFrom nlme corExp
corExp <- nlme::corExp

#' @export
#' @rdname corClasses
#' @importFrom nlme corGaus
corGaus <- nlme::corGaus

#' @export
#' @rdname corClasses
#' @importFrom nlme corLin
corLin <- nlme::corLin

#' @export
#' @rdname corClasses
#' @importFrom nlme corRatio
corRatio <- nlme::corRatio

#' @export
#' @rdname corClasses
#' @importFrom nlme corSpher
corSpher <- nlme::corSpher

#' @export
#' @rdname corClasses
#' @importFrom nlme corSymm
corSymm <- nlme::corSymm
