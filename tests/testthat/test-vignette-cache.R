test_that("vignette cache exposes cluster-level quality histograms", {
  cache_path <- file.path(find_package_root(), "inst", "extdata", "vignette_cache.rds")
  skip_if_not(file.exists(cache_path), "vignette cache is not available")

  cache <- readRDS(cache_path)

  expect_true(all(c("cache_provenance", "machine") %in% names(cache)))
  expect_true(all(c("dataset_hash", "codes_hash", "k_range", "source_hashes") %in%
                    names(cache$cache_provenance)))
  expect_identical(cache$cache_provenance$k_range, cache$k_range)
  expect_true(all(nzchar(unlist(cache$cache_provenance$source_hashes))))
  expect_false(is.null(cache$machine$inflect_version))
  expect_false(is.na(cache$machine$inflect_version))
  expect_null(cache$session_info)

  expect_true(all(c(
    "cluster_quality_df",
    "cluster_quality_summary",
    "cluster_quality_bins"
  ) %in% names(cache)))

  expected_methods <- c(
    "FlowSOM consensus",
    "FlowSOM auto-k",
    "Hierarchical (ward.D2)",
    "fastINFLECT inflection",
    "fastINFLECT threshold"
  )

  expect_true(all(expected_methods %in% cache$cluster_quality_df$method))
  expect_true(all(c("method", "k", "cluster", "cluster_unimodality") %in%
                    names(cache$cluster_quality_df)))
  expect_true(all(cache$cluster_quality_df$cluster_unimodality >= 0))
  expect_true(all(cache$cluster_quality_df$cluster_unimodality <= 100))
  expect_true(all(c("method", "quality_bin", "cluster_fraction") %in%
                    names(cache$cluster_quality_bins)))
  fractions_by_method <- tapply(
    cache$cluster_quality_bins$cluster_fraction,
    cache$cluster_quality_bins$method,
    sum
  )
  expect_equal(as.numeric(fractions_by_method), rep(1, length(expected_methods)),
               tolerance = 1e-8)
})
