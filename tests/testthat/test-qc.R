test_that("FlowSOMQC spread-only test handles constant marker values", {
  source_pkg_file("FlowSOM-QC.R")

  fake <- make_fake_flowsom(
    matrix(rep(1, 20), ncol = 2, dimnames = list(NULL, c("CD3", "CD4")))
  )

  result <- FlowSOMQC(
    FlowSOM.results = fake,
    metaclustering = as.integer(1),
    zeroes.in = TRUE,
    uniform.test = "spread",
    verbose = FALSE
  )

  expect_true(all(result, na.rm = TRUE))
  expect_equal(colnames(result), c("CD3", "CD4"))
})
