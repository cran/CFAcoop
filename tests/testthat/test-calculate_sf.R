test_that("format", {
  data("CFAdata")
  D <- subset.data.frame(
    x = CFAdata,
    subset = cell.line == levels(CFAdata$cell.line)[1]
  )
  fit.0 <- pwr_reg(seeded = D$`Cells seeded`, D$`0 Gy`)
  fit.4 <- pwr_reg(seeded = D$`Cells seeded`, D$`4 Gy`)
  expect_error(calculate_sf(
    par_ref = fit.0,
    par_treat = as.data.frame(fit.4$coefficients[, 1])
  ))
  expect_error(calculate_sf(
    par_ref = fit.0,
    par_treat = matrix(fit.4$coefficients[, 1],
      ncol = 2
    )
  ))
  expect_error(
    calculate_sf(
      par_ref = matrix(
        data = fit.0$coefficients[, 1],
        ncol = 2
      ),
      par_treat = matrix(
        data = rep(fit.4$coefficients[, 1], times = 3),
        ncol = 2
      )
    )
  )
  expect_equal(
    length(
      calculate_sf(
        par_ref = matrix(data = c(2, 1, 2, 1, 4, 3), ncol = 2, byrow = TRUE),
        par_treat = matrix(data = c(3, 2, 3, 2, 7, 6), ncol = 2, byrow = TRUE),
        c_range = c(5, 50)
      )
    ),
    3
  )
})

test_that("robust calculation", {
  data("CFAdata")
  D <- subset.data.frame(
    x = CFAdata,
    subset = cell.line == levels(CFAdata$cell.line)[1]
  )
  fit.0 <- pwr_reg(seeded = D$`Cells seeded`, D$`0 Gy`)
  fit.4 <- pwr_reg(seeded = D$`Cells seeded`, D$`4 Gy`)
  expect_equal(
    calculate_sf(par_ref = fit.0, par_treat = fit.4),
    calculate_sf(
      par_ref = matrix(
        data = fit.0$coefficients[, 1],
        ncol = 2
      ),
      par_treat = matrix(
        data = fit.4$coefficients[, 1],
        ncol = 2
      )
    )
  )
  fit.0$coefficients[2, 1] <- 0
  exp_out <- rep(Inf, 3)
  names(exp_out) <- c(5, 20, 100)
  expect_equal(calculate_sf(par_ref = fit.0, par_treat = fit.4), exp_out)
})

test_that("correct calculation", {
  data("CFAdata")
  D <- subset.data.frame(
    x = CFAdata,
    subset = cell.line == levels(CFAdata$cell.line)[1]
  )
  fit.0 <- pwr_reg(seeded = D$`Cells seeded`, D$`0 Gy`)
  fit.4 <- pwr_reg(seeded = D$`Cells seeded`, D$`4 Gy`)
  exp_out <- rep(1, 3)
  names(exp_out) <- c(5, 20, 100)
  expect_equal(
    calculate_sf(
      par_ref = matrix(data = fit.0$coefficients[, 1], ncol = 2),
      par_treat = matrix(data = fit.0$coefficients[, 1], ncol = 2)
    ),
    exp_out
  )
  expect_equal(
    calculate_sf(
      par_ref = matrix(data = c(2, 1, 2, 1), ncol = 2, byrow = TRUE),
      par_treat = matrix(data = c(3, 2, 3, 2), ncol = 2, byrow = TRUE),
      c_range = 50
    ),
    rep(exp((log(50) - 2) / 1 - (log(50) - 3) / 2), 2)
  )
})
