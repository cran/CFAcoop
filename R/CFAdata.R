#' Colony Formation Assay data on cellular cooperation
#'
#' Clonogenic survival data from seven cell lines T47D, MDA-MB231, A549,
#' HCC1806, SKBR3, SKLU1 and BT20 as presented in
#' Figure 2 in Brix et al. (2020).
#'
#' @docType data
#'
#' @usage data(CFAdata)
#'
#' @format \code{data.frame}
#'
#' @keywords dataset
#'
#' @references Brix, N., Samaga, D., Hennel, R. et al.
#' "The clonogenic assay: robustness of plating efficiency-based analysis is
#' strongly compromised by cellular cooperation." Radiat Oncol 15, 248 (2020).
#' <doi:10.1186/s13014-020-01697-y>
#'
#' @examples
#' data(CFAdata)
#' head(CFAdata)
#' cll <- levels(CFAdata$cell.line)
"CFAdata"
