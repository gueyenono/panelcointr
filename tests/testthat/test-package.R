test_that("package namespace loads", {
  expect_true("panelcointr" %in% loadedNamespaces())
})
