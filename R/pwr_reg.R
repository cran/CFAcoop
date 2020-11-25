#' @title  pwr_reg
#'
#' @description \code{pwr_reg} performs a power regression
#'   (log(C) = log(a) + b * log(S) + e)) for clonogenic assay data of
#'   experiments examining the cellular cooperation.
#'
#' @param seeded numeric vector with number of cells seeded (S)
#' @param counted numeric vector with number of colonies counted (C, same
#'   length as \code{seeded})
#'
#' @return \code{summary.lm} object as returned by \code{\link{summary}}
#'
#' @examples
#' pwr_reg(
#'   seeded = 10^(seq(1, 5, 0.5)),
#'   counted = 0.4 * (10^seq(1, 5, 0.5))^1.25 * rnorm(n = 9, 1, 0.05)
#' )
#' data(CFAdata)
#' D <- subset.data.frame(
#'   x = CFAdata,
#'   subset = cell.line == levels(CFAdata$cell.line)[1]
#' )
#' pwr_reg(seeded = D$`Cells seeded`, counted = D$`0 Gy`)
#' @export
#' @importFrom stats "lm"
pwr_reg <- function(seeded, counted) {
  if (!is.numeric(seeded) | !is.numeric(counted)) {
    stop("error: input must be numeric")
  }
  if (length(seeded) != length(counted)) {
    stop("error: input vectors must be of identical length")
  }
  x <- data.frame("S" = seeded, "C" = counted)
  x <- x[(!is.na(x$S) & !is.na(x$C)), ]
  if (sum(x <= 0) > 0) {
    warning(
      "log(0) = -Inf; power regression not applicable to null-valued variables;
      non-positive data points removed from analysis, consider regressing mean
      values."
    )
    x <- x[((x$S > 0) & (x$C > 0)), ]
  }
  if (nrow(x) < 3) {
    stop("error: not enough data for power regression")
  }
  x$lnS <- log(x$S)
  x$lnC <- log(x$C)
  fit <- lm(formula = "lnC ~ 1 + lnS", data = x)
  sfit <- summary(fit)
  return(sfit)
}
