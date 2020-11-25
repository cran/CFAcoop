test_that("input format", {
  S <- 10^seq(1, 4, 0.25)
  a <- 0.4
  b <- 1.1
  D <- data.frame(
    "S" = rep(S, each = 3),
    "C" = a * rep(S, each = 3)^b * rnorm(length(S) * 3, 1, 0.05)
  )
  expect_error(analyze_survival())
  expect_error(analyze_survival(RD = c()))
  expect_error(analyze_survival(RD = D, xtreat = 1:5))
  expect_is(analyze_survival(RD = D, xtreat = 0), "list")
})

test_that("output format", {
  S <- 10^seq(1, 4, 0.25)
  a <- 0.4
  b <- 1.1
  D <- data.frame(
    "S" = rep(S, each = 3),
    "C" = a * rep(S, each = 3)^b * rnorm(length(S) * 3, 1, 0.05)
  )
  D$C2 <- (a * 0.8) * rep(S, each = 3)^(b * 1.1) *
    rnorm(length(S) * 3, 1, 0.05)
  expect_type(analyze_survival(RD = D), "list")
  expect_length(analyze_survival(RD = D), 7)

  expect_equal(class(analyze_survival(RD = D)$name), "character")
  expect_true(is.numeric(analyze_survival(RD = D)$xtreat))
  expect_equal(class(analyze_survival(RD = D)$raw), "data.frame")
  expect_equal(class(analyze_survival(RD = D)$mean), "data.frame")
  seed.bypass <- rpois(n = 1, lambda = round(100 * rnorm(1, 20, 1)))
  seed.inner <- rpois(n = 1, lambda = 1000)
  set.seed(seed.inner)
  one <- analyze_survival(RD = D)
  set.seed(seed.inner)
  two <- analyze_survival(RD = as.matrix(D))
  expect_equal(one, two)
  set.seed(seed.bypass)
  expect_equal(class(analyze_survival(RD = D)$fit), "list")
  expect_equal(class(analyze_survival(RD = D)$fit[[1]]), "summary.lm")
  expect_equal(class(analyze_survival(RD = D)$SF), "list")
  expect_equal(class(analyze_survival(RD = D)$SF[[1]]), "numeric")
  expect_equal(class(analyze_survival(RD = D)$uncertainty), "list")
  expect_equal(class(analyze_survival(RD = D)$uncertainty[[1]])[1], "matrix")
})

test_that("find c_range", {
  S <- 10^seq(2, 4, 0.25)
  a <- 0.8
  b <- 1.1
  D <- data.frame(
    "S" = rep(S, each = 3),
    "C" = a * rep(S, each = 3)^b * rnorm(length(S) * 3, 1, 0.05)
  )
  D$C2 <- (a * 0.001) * rep(S, each = 3)^(b * 1.1) *
    rnorm(length(S) * 3, 1, 0.05)
  expect_warning(analyze_survival(RD = D))
})
