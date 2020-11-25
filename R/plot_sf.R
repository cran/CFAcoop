#' @title  plot_sf
#'
#' @description plot cellular cooperativity and clonogenic survival for
#' colony formation assay data
#'
#' @param SF list build of objects returned by \code{analyze_survival}
#' @param showUncertainty logical, switches on/off uncertainty bands for
#' sf-values.
#' @return none
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
#' plot_sf(SF)
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
#' plot_sf(SF)
#' @importFrom grDevices "col2rgb" "colorRampPalette" "rgb"
#' @importFrom graphics "abline" "axis" "par" "plot" "title"
#' @importFrom Hmisc "errbar"
#' @importFrom graphics "polygon"
#' @export
#'
plot_sf <- function(SF, showUncertainty = TRUE) {
  if (length(SF) > 10) {
    stop(
      "error: more than ten experiments were chosen for plotting.
      Consider separating data set for presentation."
    )
  } else if (length(SF) == 0) {
    stop(
      "error: empty SF object cannot be plotted."
    )
  }
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(2, length(SF)))
  par(
    mar = c(2.5, 3.5, 0.5, 0.5),
    mgp = c(1.5, 0.5, 0)
  )
  collect_sf <- data.frame(
    "Exp" = NA,
    "treat" = NA,
    "sf" = NA,
    "q1" = NA,
    "q2" = NA
  )
  alpha_bg <- 2 * 42
  for (t in seq_along(SF)) {
    CurSF <- SF[[t]]
    N_treat <- length(CurSF$fit)
    if (CurSF$name == "no name") {
      CurSF$name <- paste0("Experiment ", t)
    }
    CurSF$plot <- CurSF$raw
    CurSF$plot[CurSF$plot == 0] <- 0.1
    CurSF$plot <- log(CurSF$plot) / log(10)
    # log_b(x) = log_a(x)/log_a(b) for plotting
    CurSF$pm <- log(CurSF$mean) / log(10)
    CurSF$pm[CurSF$pm == -Inf] <- NaN

    with_data <-
      apply(
        X = !is.na(CurSF$raw[, c(2:(dim(CurSF$raw)[2]))]),
        MARGIN = 1,
        FUN = "sum",
        na.rm = TRUE
      ) >= 1
    x_lim <- range(CurSF$plot[with_data, 1])
    y_lim <- c(-1, 3)
    colhex <- colorRampPalette(c("#43E08700", "#00612A00"))(N_treat) # D23264
    colors <- col2rgb(colhex) / 255
    alpha <- 0.42
    plot(
      x = CurSF$plot[, 1],
      y = CurSF$plot[, 2],
      main = "",
      xlim = x_lim,
      ylim = y_lim,
      pch = "+",
      yaxt = "n",
      xlab = "cells seeded",
      ylab = "colonies counted",
      axes = FALSE,
      col = rgb(
        red = colors[1, 1],
        green = colors[2, 1],
        blue = colors[3, 1],
        alpha = alpha,
        maxColorValue = 1
      )
    )
    xtick <- floor(min(CurSF$plot[, 1], na.rm = TRUE)):
    ceiling(max(CurSF$plot[, 1], na.rm = TRUE))
    ytick <- -1:3
    axis(
      side = 2,
      at = ytick,
      las = 1,
      labels = c("no cols.", "1", "10", "100", "1000")
    )
    axis(
      side = 1,
      at = xtick,
      las = 1,
      labels = 10^xtick
    )
    par(new = TRUE)
    polygon(
      x = c(x_lim * c(0.8, 1.2), rev(x_lim * c(0.8, 1.2))),
      y = log10(c(5, 5, 100, 100)),
      border = NA,
      col = rgb(25, 25, 25, alpha = 42, maxColorValue = 255)
    )
    par(new = TRUE)
    plot(
      x = CurSF$pm[, 1],
      y = CurSF$pm[, 3],
      xlim = x_lim,
      ylim = y_lim,
      pch = 19,
      ann = FALSE,
      xaxt = "n",
      yaxt = "n",
      axes = FALSE,
      col = rgb(
        red = colors[1, 1],
        green = colors[2, 1],
        blue = colors[3, 1],
        alpha = 1,
        maxColorValue = 1
      )
    )

    abline(
      a = 0, b = 1, lty = 2,
      col = rgb(25, 25, 25, alpha = alpha_bg, maxColorValue = 255)
    )
    abline(
      a = -1, b = 1, lty = 2,
      col = rgb(25, 25, 25, alpha = alpha_bg, maxColorValue = 255)
    )
    abline(
      a = -2, b = 1, lty = 2,
      col = rgb(25, 25, 25, alpha = alpha_bg, maxColorValue = 255)
    )
    abline(
      a = -3, b = 1, lty = 2,
      col = rgb(25, 25, 25, alpha = alpha_bg, maxColorValue = 255)
    )
    abline(
      a = -4, b = 1, lty = 2,
      col = rgb(25, 25, 25, alpha = alpha_bg, maxColorValue = 255)
    )
    abline(
      a = CurSF$fit[[1]]$coefficients[1, 1] / log(10),
      b = CurSF$fit[[1]]$coefficients[2, 1],
      col = rgb(
        red = colors[1, 1],
        green = colors[2, 1],
        blue = colors[3, 1],
        alpha = 1,
        maxColorValue = 1
      )
    )
    collect_sf <- rbind(collect_sf, c(t, CurSF$xtreat[1], 1, 1, 1))
    sf_vec <- NULL
    q1_vec <- NULL
    q2_vec <- NULL
    x_vec <- NULL
    col_vec <- colhex[1]
    for (i in 2:N_treat) {
      par(new = TRUE)
      plot(
        x = CurSF$plot[, 1],
        y = CurSF$plot[, i + 1],
        xlim = x_lim,
        ylim = y_lim,
        pch = "+",
        ann = FALSE,
        axes = FALSE,
        xaxt = "n",
        yaxt = "n",
        col = rgb(
          red = colors[1, i],
          green = colors[2, i],
          blue = colors[3, i],
          alpha = alpha,
          maxColorValue = 1
        )
      )
      par(new = TRUE)
      plot(
        x = CurSF$pm[, 1],
        y = CurSF$pm[, i + 2],
        xlim = x_lim,
        ylim = y_lim,
        pch = 19,
        ann = FALSE,
        xaxt = "n",
        yaxt = "n",
        axes = FALSE,
        col = rgb(
          red = colors[1, i],
          green = colors[2, i],
          blue = colors[3, i],
          alpha = 1,
          maxColorValue = 1
        )
      )
      abline(
        a = CurSF$fit[[i]]$coefficients[1, 1] / log(10),
        b = CurSF$fit[[i]]$coefficients[2, 1],
        col = rgb(
          red = colors[1, i],
          green = colors[2, i],
          blue = colors[3, i],
          alpha = 1,
          maxColorValue = 1
        )
      )
      sf_vec <- c(sf_vec, CurSF$SF[[i - 1]])
      q1_vec <- c(q1_vec, CurSF$uncertainty[[i - 1]][, 1])
      q2_vec <- c(q2_vec, CurSF$uncertainty[[i - 1]][, 2])
      x_vec <- c(x_vec, rep(CurSF$xtreat[i], length(CurSF$SF[[i - 1]])))
      col_vec <- c(col_vec, rep(colhex[i], length(CurSF$SF[[i - 1]])))
    }
    keep_sf <-
      data.frame(
        "Exp" = rep(t, length(x_vec)),
        "treat" = x_vec,
        "sf" = sf_vec,
        "q1" = q1_vec,
        "q2" = q2_vec
      )
    collect_sf <- rbind(collect_sf, keep_sf)
  }
  collect_sf <- collect_sf[-1, ]

  par(mar = c(2.5, 3.75, 2.5, 0.5), mgp = c(1.5, 0.5, 0))
  for (sfi in seq_along(SF)) {
    if (SF[[sfi]]$name == "no name") {
      SF[[sfi]]$name <- paste0("Experiment ", sfi)
    }
    PD <-
      subset.data.frame(x = collect_sf, subset = (collect_sf$"Exp" == sfi))
    plot(
      x = PD$treat,
      log10(PD$sf),
      main = SF[[sfi]]$name,
      col.main = rgb(
        red = 0, green = 148, blue = 64, alpha = 255,
        maxColorValue = 255
      ),
      las = 1,
      ylab = "",
      xaxt = "n",
      yaxt = "n",
      xlab = "treatment",
      ylim = range(log10(collect_sf$sf)),
      col = col_vec,
      axes = FALSE,
      pch = 19
    )
    ytick <- round(min(log10(collect_sf$sf), na.rm = TRUE)):
    round(max(log10(collect_sf$sf), na.rm = TRUE))
    axis(
      side = 2,
      at = ytick,
      las = 1,
      labels = paste0(10^ytick * 100, "%")
    )
    title(
      ylab = "clonogenic survival",
      line = 2.5
    )
    axis(
      side = 1,
      at = PD$treat,
      labels = PD$treat
    )
    if (showUncertainty) {
      with(data = PD, errbar(treat, log10(sf), log10(q1), log10(q2),
        col = col_vec,
        add = TRUE, pch = 1, errbar.col = col_vec
      ))
    }
  }
}
