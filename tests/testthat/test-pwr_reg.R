test_that("reject wrong input", {
  expect_error(pwr_reg(seeded = letters[1:5],counted = 1:5))
  expect_error(pwr_reg(seeded = c(),counted = rnorm(5)))
  expect_error(pwr_reg(seeded = 15,counted=5))
  expect_error(pwr_reg(seeded = 1:15,counted=1:10))
})

test_that("dealing with NA", {
  noise <- rnorm(11,1,0.05)
  expect_equal(pwr_reg(seeded = c(1:5,NA,6:10),
                           counted = c(1:5,5:10)*noise),
                   pwr_reg(seeded = c(1:5,5:10),
                           counted = c(1:5,NA,6:10)*noise))
})

test_that("power regression boundaries", {
  expect_warning(pwr_reg(seeded = -1:5,counted = rnorm(7,mean = 5)))
  S <- 10^seq(1,4,1/3)
  C <- 0.42 * S^1.2
  expect_warning(pwr_reg(seeded = S,counted = C))
})

test_that("output format", {
  expect_identical(class(pwr_reg(seeded = 1:10,counted = rnorm(n = 10,10))),
                   "summary.lm")
  expect_type(pwr_reg(seeded = 1:10,counted = rnorm(n = 10,10)),"list")
  expect_length(pwr_reg(seeded = 1:10,counted = rnorm(n = 10,10)),11)
})

test_that("unbiased regression", {
  S <- 10^seq(1,4,1/3)
  a <- 0.42
  b <- 1.2
  speicher <- rep(NA,100)
  for (i in seq_along(speicher)){
    C <- a * S^b * rnorm(n = length(S),mean = 1,sd = 0.05)
    speicher[i] <- pwr_reg(seeded = S,counted = C)$coefficients[2,1]
  }
  expect_gt(object = mean(speicher),expected = b - 5*sd(speicher)/10)
  expect_lt(object = mean(speicher),expected = b + 5*sd(speicher)/10)
})
