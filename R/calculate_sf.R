#' @title  calculate_sf
#'
#' @description calculates the survival fraction according
#'   to the procedure presented in Brix et al. (2020), which is robust against
#'   cellular cooperation.
#'
#' @param par_ref \code{summary.lm} object or 2-column matrix for the
#'   treatment-free reference survival
#' @param par_treat \code{summary.lm} object or 2-column matrix for the
#'   clonogenic survival after treatment
#' @param C colony number for which the survival fraction is calculated
#'   (default = 20)
#'
#' @return survival fractions.
#'   If par_ref and par_treat are \code{summary.lm} objects,
#'   a scalar is returned.
#'   If par_ref and par_treat are matrices,
#'   a vector of the same length as nrow(par_treat) is returned
#'
#' @examples
#' seeded <- 10^(seq(1, 5, 0.5))
#' counted.ref <- 0.4 * 10^(seq(1, 5, 0.5) + rnorm(n = 9, 0, 0.1))^1.1
#' counted.treat <- 0.01 * 10^(seq(1, 5, 0.5) + rnorm(n = 9, 0, 0.1))^1.2
#' fit_ref <- pwr_reg(seeded = seeded, counted = counted.ref)
#' fit_treat <- pwr_reg(seeded = seeded, counted = counted.treat)
#' calculate_sf(par_ref = fit_ref, par_treat = fit_treat)
#' data("CFAdata")
#' D <- subset.data.frame(
#'   x = CFAdata,
#'   subset = cell.line == levels(CFAdata$cell.line)[1]
#' )
#' fit_ref <- pwr_reg(seeded = D$`Cells seeded`, counted = D$`0 Gy`)
#' fit_treat <- pwr_reg(seeded = D$`Cells seeded`, counted = D$`4 Gy`)
#' calculate_sf(par_ref = fit_ref, par_treat = fit_treat)
#' @export
#'
calculate_sf <- function(par_ref, par_treat, C = 20) {
  if (!prod(c(class(par_ref)[1], class(par_treat)[1]) %in% c(
    "summary.lm",
    "matrix"
  ))) {
    stop("error: par_ref and par_treat must be of class summary.lm or matrix")
  }
  if (class(par_ref)[1] != class(par_treat)[1]) {
    stop("error: class of par_ref and par_treat must be identical ")
  }
  if (class(par_ref)[1] == "summary.lm") {
    # calculate survival fraction from two sfit-objects or two pairs c(a,b)
    SF <-
      exp(((log(C) - par_ref$coefficients[1, 1]) /
        par_ref$coefficients[2, 1]) -
        ((log(C) - par_treat$coefficients[1, 1]) /
          par_treat$coefficients[2, 1]))
    names(SF) <- C
  } else {
    if (!identical(dim(par_ref), dim(par_treat))) {
      stop("error: par_ref and par_treat must be of identical size")
    }
    if (nrow(par_ref) > 1) {
      SF <- exp(((log(C) - par_ref[, 1]) / par_ref[, 2]) -
        ((log(C[1]) - par_treat[, 1]) / par_treat[, 2]))
    } else {
      SF <- exp(((log(C) - par_ref[1]) / par_ref[2]) -
        ((log(C) - par_treat[1]) / par_treat[2]))
      names(SF) <- C
    }
  }
  return(SF)
}
