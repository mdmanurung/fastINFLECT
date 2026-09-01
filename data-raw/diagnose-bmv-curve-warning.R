#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) {
  stop("Could not resolve the diagnostic script path.", call. = FALSE)
}
script_file <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
evidence_dir <- file.path(
  repo_root,
  "inst",
  "benchmarks",
  "real-model-selection"
)
selection_path <- file.path(evidence_dir, "selection-evidence.rds")
if (!file.exists(selection_path)) {
  stop("Selection evidence is incomplete: ", selection_path, call. = FALSE)
}

atomic_save_rds <- function(object, path) {
  temporary <- tempfile("curve-warning-", tmpdir = dirname(path))
  saveRDS(object, temporary, compress = TRUE)
  if (!file.rename(temporary, path)) {
    stop("Could not publish ", path, call. = FALSE)
  }
  invisible(path)
}

atomic_write_csv <- function(object, path) {
  temporary <- tempfile("curve-warning-", tmpdir = dirname(path))
  utils::write.csv(object, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) {
    stop("Could not publish ", path, call. = FALSE)
  }
  invisible(path)
}

suppressWarnings(pkgload::load_all(repo_root, quiet = TRUE))
selection <- readRDS(selection_path)
summary <- selection$full_inflect$criterion_summary
stored <- selection$full_inflect$criterion_selections
criteria <- unique(summary$criterion)

rows <- lapply(criteria, function(criterion) {
  data <- summary[
    summary$criterion == criterion,
    c("k", "qc_pass_rate", "criterion"),
    drop = FALSE
  ]
  warning_messages <- character()
  warning_calls <- character()
  fit <- withCallingHandlers(
    fastINFLECT::QC.to.curve(data),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      warning_calls <<- c(
        warning_calls,
        paste(deparse(conditionCall(w)), collapse = " ")
      )
      invokeRestart("muffleWarning")
    }
  )
  expected_k <- stored$k[
    stored$criterion == criterion &
      stored$method == "inflection"
  ]
  data.frame(
    criterion = criterion,
    stored_inflection_k = expected_k,
    rerun_inflection_k = as.integer(fit$lfunction$knee),
    warning_count = length(warning_messages),
    warning_messages = paste(unique(warning_messages), collapse = " | "),
    warning_calls = paste(unique(warning_calls), collapse = " | "),
    all_fitted_values_finite = all(is.finite(
      fit$fittedcurve$qc_pass_rate
    )),
    observed_range_fitted_values_finite = all(is.finite(
      fit$fittedcurve$qc_pass_rate[
        fit$fittedcurve$k %in% data$k
      ]
    )),
    inflection_matches = identical(
      as.integer(fit$lfunction$knee),
      as.integer(expected_k)
    ),
    stringsAsFactors = FALSE
  )
})
diagnostic <- do.call(rbind, rows)
rownames(diagnostic) <- NULL

expected_call <- grepl(
  "log[(]dose/parmMat",
  diagnostic$warning_calls,
  fixed = FALSE
)
classification <- if (
  all(diagnostic$warning_messages == "NaNs produced") &&
    all(expected_call) &&
    all(diagnostic$all_fitted_values_finite) &&
    all(diagnostic$observed_range_fitted_values_finite) &&
    all(diagnostic$inflection_matches)
) {
  "transient_LL4_optimizer_exploration_warning"
} else {
  "unresolved_curve_fit_warning"
}
if (!identical(classification, "transient_LL4_optimizer_exploration_warning")) {
  stop("Curve-fit warning diagnostic did not pass.", call. = FALSE)
}

bundle <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%OS%z"),
  classification = classification,
  interpretation = paste0(
    "drc::drm explored a temporary nonpositive LL.4 ED50 and warned while ",
    "evaluating log(dose / e). Final fitted values were finite and all ",
    "criterion inflections exactly matched the stored run. The warning is ",
    "retained in the selection runtime evidence."
  ),
  diagnostic = diagnostic,
  selection_file = normalizePath(selection_path, mustWork = TRUE),
  selection_md5 = unname(tools::md5sum(selection_path)),
  script_file = script_file,
  script_md5 = unname(tools::md5sum(script_file)),
  dependency_versions = c(
    R = as.character(getRversion()),
    fastINFLECT = as.character(utils::packageVersion("fastINFLECT")),
    drc = as.character(utils::packageVersion("drc"))
  ),
  session_info = utils::sessionInfo()
)
atomic_save_rds(
  bundle,
  file.path(evidence_dir, "curve-fit-warning-diagnostic.rds")
)
atomic_write_csv(
  diagnostic,
  file.path(evidence_dir, "curve-fit-warning-diagnostic.csv")
)
message("Curve-fit warning diagnostic passed: ", classification)
