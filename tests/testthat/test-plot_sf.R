test_that("input", {
  data("CFAdata")
  SF.list <- vector("list", 14)
  for (i in 1:7) {
    SF.list[[i]] <- analyze_survival(
      RD = subset.data.frame(
        x = CFAdata,
        subset = (CFAdata$cell.line == levels(CFAdata$cell.line)[i])
      )[, -1]
    )
    SF.list[[7 + i]] <- analyze_survival(
      RD = subset.data.frame(
        x = CFAdata,
        subset = (CFAdata$cell.line == levels(CFAdata$cell.line)[i])
      )[, -1]
    )
  }
  expect_error(plot_sf(SF = SF.list))
  expect_error(plot_sf(SF = c()))
})

test_that("no output", {
  data("CFAdata")
  SF.list <- vector("list", 7)
  for (i in seq_along(SF.list)) {
    SF.list[[i]] <- analyze_survival(
      RD = subset.data.frame(
        x = CFAdata,
        subset = (CFAdata$cell.line == levels(CFAdata$cell.line)[i])
      )[, -1]
    )
  }
  expect_equal(class(plot_sf(SF = SF.list)), "NULL")
})
