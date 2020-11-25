#' @title  analyze_survival
#'
#' @description wrapper function for robust analysis clonogenic survival data
#'   from the colony formation assay according to Brix et al. (2020),
#'   Radiation Oncology.
#'   Mean values are calculated and used for power regression.
#'   Resulting coefficients are used for
#'   calculation of survival fractions and corresponding uncertainty analysis.
#'
#' @param RD data.frame or matrix containing a table of experiment data
#' @param name optional: experiment name (e.g. name of cell line)
#' @param xtreat optional: treatment dose of the colonies counted in the
#'   corresponding columns of RD
#' @param c_range number or vector of numbers of colonies counted for which
#'   the survival fraction is to be calculated (default = c(5, 20, 100))
#'
#' @return list object containing several experiments and treatments organized
#'   for convenient plotting with \code{plot_sf}
#'
#' @examples
#' seeded <- rep(10^(seq(1,5,0.5)),each = 3)
#' df.1 <- data.frame(
#'   "seeded" = seeded,
#'   "counted1" = 0.4 * seeded^1.1 * rnorm(n = length(seeded),1,0.05),
#'   "counted2" = 0.2 * seeded^1.125 * rnorm(n = length(seeded),1,0.05),
#'   "counted3" = 0.05 * seeded^1.25 * rnorm(n = length(seeded),1,0.05))
#' df.2 <- data.frame("seeded" = seeded,
#'   "counted1" = 0.5 * seeded^1.01 * rnorm(n = length(seeded),1,0.05),
#'   "counted2" = 0.4 * seeded^1.0125 * rnorm(n = length(seeded),1,0.05),
#'   "counted3" = 0.2 * seeded^1.025 * rnorm(n = length(seeded),1,0.05))
#' SF <- vector("list",2)
#' SF[[1]] <- analyze_survival(RD = df.1,
#'                             name = "cell line a",
#'                             xtreat = c(0,1,4),
#'                             c_range = c(5,20,100))
#' SF[[2]] <- analyze_survival(RD = df.2,
#'                             name = "cell line b",
#'                             xtreat = c(0,1,4))
#' @importFrom stats "aggregate" "quantile" "vcov"
#' @importFrom mvtnorm "rmvnorm"
#' @export
#'
analyze_survival <- function(RD,
                             name = "no name",
                             xtreat = NULL,
                             c_range = c(5,20,100)) {
  if (!(class(RD)[1] %in% c("data.frame","matrix"))){
    stop("error: RD must be of class data.frame or matrix")
  }
  result <- list("name" = name)
  if (is.null(xtreat)) {
    result$"xtreat" <- 0:(ncol(RD) - 2)
  } else {
    if (length(xtreat) != (ncol(RD) - 1)) {
      stop("error: length of assigned treatments does not match data ")
    }
    result$"xtreat" <- xtreat
  }
  result$"raw" <- as.data.frame(RD)
  result$"mean" <-
    aggregate(
      x = result$raw,
      by = list(result$raw[, 1]),
      FUN = "mean",
      na.rm = TRUE
    )
  result$"fit" <- vector("list", dim(RD)[2] - 1)
  result$"SF" <- vector("list", dim(RD)[2] - 2)
  result$"uncertainty" <- vector("list", dim(RD)[2] - 2)
  result$"fit"[[1]] <- pwr_reg(seeded = result$"mean"[, 2],
                               counted = result$"mean"[, 3])
  for (i in seq_along(result$"fit")[-1]) {
    result$"fit"[[i]] <- pwr_reg(seeded = result$"mean"[, 2],
                                 counted = result$"mean"[, i + 2])
    ran_ref <- range(result$"raw"[, 2],na.rm = TRUE)
    ran_tre <- range(result$"raw"[, i + 1],na.rm = TRUE)
    if ( (ran_ref[2] >= ran_tre[1] ) & (ran_ref[1]) <= ran_tre[2] ){
      result$"SF"[[i - 1]] <- calculate_sf(
        par_ref = result$"fit"[[1]],
        par_treat = result$"fit"[[i]],
        c_range = c_range
      )
      Pref<- mvtnorm::rmvnorm(n = 10^3,
                     mean = result$fit[[1]]$coefficients[,1],
                     sigma = vcov(result$fit[[1]]))
      Ptreat<- mvtnorm::rmvnorm(n = 10^3,
                       mean = result$fit[[i]]$coefficients[,1],
                       sigma = vcov(result$fit[[i]]))
      result$"uncertainty"[[i - 1]] <- matrix(
        data = rep(NA,2*length(c_range)),
        ncol = 2)
      for (j in seq_along(c_range)){
        result$"uncertainty"[[i - 1]][j,] <- quantile(
          x = calculate_sf(par_ref = Pref,
                           par_treat = Ptreat,
                           c_range = c_range[j]),
          probs = c(0.025,0.975))
      }
      rownames(result$"uncertainty"[[i - 1]]) <-
        names(result$"SF"[[i - 1]])
    } else {
      warning("warning: SF calculation omitted, range of colonies counted
              in reference and treated cells do not overlap.")
      result$"SF"[[i - 1]] <- NaN
      result$"uncertainty"[[i - 1]] <- NaN
    }
  }
  return(result)
}
