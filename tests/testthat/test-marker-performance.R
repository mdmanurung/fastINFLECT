test_that("marker.performance works with omitted and single marker selections", {
  source_pkg_file("plot-markerperformance.R")

  accuracy <- list(
    "5" = matrix(
      c(TRUE, FALSE, TRUE, TRUE),
      nrow = 2,
      dimnames = list(c("cluster1", "cluster2"), c("CD3", "CD4"))
    ),
    "6" = matrix(
      c(TRUE, TRUE, FALSE, TRUE),
      nrow = 2,
      dimnames = list(c("cluster1", "cluster2"), c("CD3", "CD4"))
    )
  )
  result <- structure(list(accuracy.sets = accuracy), class = "inflect.results")

  all_markers <- marker.performance(result)
  single_marker <- marker.performance(result, markers = "CD3")

  expect_s3_class(all_markers$plot, "ggplot")
  expect_equal(sort(unique(as.character(all_markers$marker.dataframe$Marker))), c("CD3", "CD4"))
  expect_equal(unique(as.character(single_marker$marker.dataframe$Marker)), "CD3")
})
