test_that("benchmark cache exposes legacy timing and marker-expression panels", {
  cache_path <- file.path(find_package_root(), "inst", "extdata", "benchmark_cache.rds")
  skip_if_not(file.exists(cache_path), "benchmark cache is not available")

  cache <- readRDS(cache_path)

  expect_true(all(c("legacy_inflect", "marker_histogram", "marker_histogram_summary") %in% names(cache)))
  expect_true(all(c("legacy_inflect_scan_seconds", "speedup_legacy_inflect") %in% names(cache$totals)))
  expect_true(all(c("cache_provenance", "machine") %in% names(cache)))
  expect_true(all(c("dataset_hash", "codes_hash", "k_range", "source_hashes") %in%
                    names(cache$cache_provenance)))
  expect_identical(cache$cache_provenance$k_range, cache$k_range)
  expect_true(all(nzchar(unlist(cache$cache_provenance$source_hashes))))
  expect_false(is.null(cache$machine$inflect_version))
  expect_false(is.na(cache$machine$inflect_version))

  expect_true(is.data.frame(cache$legacy_inflect))
  expect_true(all(c("k", "seconds", "seconds_se", "unimodality", "method") %in%
                    names(cache$legacy_inflect)))
  expect_true(all(is.finite(cache$legacy_inflect$seconds)))
  expect_true(all(is.finite(cache$legacy_inflect$seconds_se)))
  expect_true(cache$totals$legacy_inflect_scan_seconds > cache$totals$inflect_scan_seconds)
  expect_true(cache$totals$legacy_inflect_scan_seconds_se > 0)
  expect_true(cache$totals$speedup_legacy_inflect > 1)
  expect_equal(cache$efficiency$legacy_projection$timing_type,
               "estimated_from_sampled_diptest_cells")
  expect_true(cache$efficiency$legacy_projection$sampled_tests <
                cache$efficiency$total_cluster_marker_tests_legacy)
  expect_equal(cache$totals$consensus_scan_seconds, sum(cache$consensus_df$seconds))
  expect_true(all(cache$consensus_df$timing_type == "partition_plus_qc_scoring"))
  expect_true(all(cache$amortized_df$timing_type == "single_consensus_run_plus_qc_scoring"))

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
