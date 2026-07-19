test_that("benchmark cache exposes legacy timing and marker-expression panels", {
  cache_path <- file.path(find_package_root(), "inst", "extdata", "benchmark_cache.rds")
  skip_if_not(file.exists(cache_path), "benchmark cache is not available")

  cache <- readRDS(cache_path)

  expect_true(all(c("legacy_inflect", "marker_histogram", "marker_histogram_summary") %in% names(cache)))
  expect_true(all(c("legacy_inflect_scan_seconds", "speedup_legacy_inflect") %in% names(cache$totals)))

  expect_true(is.data.frame(cache$legacy_inflect))
  expect_true(all(c("k", "seconds", "unimodality", "method") %in% names(cache$legacy_inflect)))
  expect_true(all(is.finite(cache$legacy_inflect$seconds)))
  expect_true(cache$totals$legacy_inflect_scan_seconds > cache$totals$inflect_scan_seconds)
  expect_true(cache$totals$speedup_legacy_inflect > 1)

  expected_histogram_cols <- c(
    "method", "k", "cluster", "marker", "expression", "qc_pass",
    "dip_pvalue", "iqr", "n_events"
  )
  expect_true(all(expected_histogram_cols %in% names(cache$marker_histogram)))
  expect_true(all(c("method", "k", "cluster", "marker", "qc_pass", "n_events") %in%
                    names(cache$marker_histogram_summary)))
  expect_true(all(cache$marker_histogram$n_events >= 5))
  expect_true(any(cache$marker_histogram$method == "fastINFLECT threshold"))
  expect_true(any(cache$marker_histogram$method == "FlowSOM consensus"))
  expect_true(any(cache$marker_histogram$qc_pass))
  expect_true(any(!cache$marker_histogram$qc_pass))
})
