#' @title  export_sf
#'
#' @description export table with results of clonogenic survival analysis
#' from the colony formation assay considering cellular cooperation
#'
#' @param SF list build of objects returned by \code{analyze_survival}
#' @return data.frame containing all estimated coefficients and effects from
#' all experiments contained in \code{SF}
#'
#' @examples
#' seeded <- rep(10^(seq(1, 5, 0.5)), each = 3)
#' df.1 <- data.frame(
#'   "seeded" = seeded,
#'   "counted1" = 0.4 * seeded^1.1 * rnorm(n = length(seeded), 1, 0.05),
#'   "counted2" = 0.2 * seeded^1.125 * rnorm(n = length(seeded), 1, 0.05),
#'   "counted3" = 0.05 * seeded^1.25 * rnorm(n = length(seeded), 1, 0.05)
#' )
#' df.2 <- data.frame(
#'   "seeded" = seeded,
#'   "counted1" = 0.5 * seeded^1.01 * rnorm(n = length(seeded), 1, 0.05),
#'   "counted2" = 0.4 * seeded^1.0125 * rnorm(n = length(seeded), 1, 0.05),
#'   "counted3" = 0.2 * seeded^1.025 * rnorm(n = length(seeded), 1, 0.05)
#' )
#' SF <- vector("list", 2)
#' SF[[1]] <- analyze_survival(
#'   RD = df.1, name = "cell line a",
#'   xtreat = c(0, 1, 4)
#' )
#' SF[[2]] <- analyze_survival(
#'   RD = df.2, name = "cell line b",
#'   xtreat = c(0, 1, 4)
#' )
#' export_sf(SF)
#'
#' data("CFAdata")
#' SF <- vector("list", 4)
#' ll <- levels(CFAdata$cell.line)[c(1, 3, 5, 7)]
#' for (i in seq_along(ll)) {
#'   cdat <- subset.data.frame(
#'     x = CFAdata,
#'     subset = CFAdata$cell.line == ll[i]
#'   )
#'   SF[[i]] <- analyze_survival(
#'     RD = cdat[, -1],
#'     name = ll[i],
#'     xtreat = c(0, 1, 2, 4, 6, 8)
#'   )
#' }
#' export_sf(SF)
#' @export
#'
export_sf <- function(SF) {
  if (class(SF) != "list") {
    stop("error: SF must be of class 'list'")
  }
  if (!is.list(SF[[1]])) {
    if (!identical(
      x = names(SF),
      y = c(
        "name", "xtreat", "raw", "mean", "fit", "SF", "uncertainty"
      )
    )) {
      stop("error: SF object is not of the form as returned by
           analyze_survival.")
    }
    SFinput <- SF
    SF <- vector("list", 1)
    SF[[1]] <- SFinput
  }
  result <- data.frame(
    "cell.line" = NA,
    "treatment" = NA,
    "ln(a)" = NA,
    "ln(a).se" = NA,
    "b" = NA,
    "b.se" = NA,
    "r" = NA
  )
  cur.index <- 1
  for (i in seq_along(SF)) {
    cSF <- SF[[i]]
    for (j in seq_along(cSF$xtreat)) {
      result[cur.index, 1] <- as.character(cSF$name[1])
      cur.exp <- c(
        cSF$xtreat[j],
        round(cSF$fit[[j]]$coefficients[1, 1], digits = 3),
        round(cSF$fit[[j]]$coefficients[1, 2], digits = 3),
        round(cSF$fit[[j]]$coefficients[2, 1], digits = 3),
        round(cSF$fit[[j]]$coefficients[2, 2], digits = 3),
        round(vcov(cSF$fit[[j]])[1, 2] /
          prod(sqrt(diag(vcov(cSF$fit[[j]])))), digits = 3)
      )
      result[cur.index, 2:7] <- cur.exp
      if (j > 1) {
        for (k in seq_along(cSF$SF[[j - 1]])) {
          SF_C <- paste0("SF(", names(cSF$SF[[j - 1]])[k], ")")
          SF_C_lb <- paste0("SF(", names(cSF$SF[[j - 1]])[k], ")_uncert_lb")
          SF_C_ub <- paste0("SF(", names(cSF$SF[[j - 1]])[k], ")_uncert_ub")
          if (!(SF_C %in% colnames(result))) {
            sofar <- colnames(result)
            result <- cbind(result, NA)
            result <- cbind(result, NA)
            result <- cbind(result, NA)
            colnames(result) <- c(sofar, c(SF_C, SF_C_lb, SF_C_ub))
          }
          result[cur.index, colnames(result) %in% SF_C] <- round(
            cSF$SF[[j - 1]][k],
            digits = 4
          )
          result[cur.index, colnames(result) %in% SF_C_lb] <- round(
            cSF$uncertainty[[j - 1]][k, 1],
            digits = 4
          )
          result[cur.index, colnames(result) %in% SF_C_ub] <- round(
            cSF$uncertainty[[j - 1]][k, 2],
            digits = 4
          )
        }
      }
      cur.index <- cur.index + 1
    }
  }
  return(result)
}
