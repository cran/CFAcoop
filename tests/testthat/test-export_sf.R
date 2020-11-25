test_that("input", {
  data("CFAdata")
  SF.list <- vector("list", 7)
  for (i in 1:7) {
    SF.list[[i]] <- analyze_survival(
      RD = subset.data.frame(
        x = CFAdata,
        subset = (CFAdata$cell.line == levels(CFAdata$cell.line)[i])
      )[, -1]
    )
  }
  expect_error(export_sf(SF = c()))
  expect_error(export_sf(SF = 42))
  SF.test <- SF.list[[1]]
  expect_equal(
    object = class(export_sf(SF = SF.test)),
    expected = "data.frame"
  )
  SF.test$SF <- NULL
  expect_error(
    object = class(export_sf(SF = SF.test)),
    expected = "data.frame"
  )
})

test_that("ouput", {
  data("CFAdata")
  SF.list <- vector("list", 7)
  for (i in 1:7) {
    SF.list[[i]] <- analyze_survival(
      RD = subset.data.frame(
        x = CFAdata,
        subset = (CFAdata$cell.line == levels(CFAdata$cell.line)[i])
      )[, -1]
    )
  }
  expect_equal(
    object = class(export_sf(SF = SF.list)),
    expected = "data.frame"
  )
})
