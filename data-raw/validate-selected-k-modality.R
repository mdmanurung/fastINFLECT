#!/usr/bin/env Rscript

## Supplemental selected-k audit for the EXV BMV model.
##
## The primary modality workflow was frozen around candidates inferred from a
## historical spread-only sweep. The dense complete-distribution run selected
## k = 27/44 for the combined criterion, k = 28/44 for dip, and k = 38/55 for
## IQR spread. k = 55 is already present in the primary audit; this wrapper
## reuses its tested implementation to add k = 27, 28, 38, and 44 without
## changing the source of an in-flight primary workflow.

wrapper_args <- commandArgs(trailingOnly = FALSE)
wrapper_file <- sub(
  "^--file=",
  "",
  wrapper_args[startsWith(wrapper_args, "--file=")]
)
if (length(wrapper_file) != 1L) {
  stop("Could not resolve the selected-k wrapper path.", call. = FALSE)
}
wrapper_file <- normalizePath(wrapper_file[[1L]], mustWork = TRUE)
wrapper_root <- normalizePath(
  file.path(dirname(wrapper_file), ".."),
  mustWork = TRUE
)
core_file <- normalizePath(
  file.path(wrapper_root, "data-raw", "validate-real-model-modality.R"),
  mustWork = TRUE
)

if (!nzchar(Sys.getenv("INFLECT_MODALITY_WORKDIR"))) {
  Sys.setenv(
    INFLECT_MODALITY_WORKDIR = file.path(
      wrapper_root,
      "data-raw",
      ".selected-k-modality-work"
    )
  )
}
if (!nzchar(Sys.getenv("INFLECT_MODALITY_EVIDENCE_DIR"))) {
  Sys.setenv(
    INFLECT_MODALITY_EVIDENCE_DIR = file.path(
      wrapper_root,
      "inst",
      "benchmarks",
      "selected-k-modality"
    )
  )
}

core_lines <- readLines(core_file, warn = FALSE)
switch_index <- grep("^switch\\($", core_lines)
if (length(switch_index) != 1L) {
  stop("Could not isolate the modality workflow entry point.", call. = FALSE)
}
core_code <- paste(
  core_lines[seq_len(switch_index[[1L]] - 1L)],
  collapse = "\n"
)
core <- new.env(parent = globalenv())
eval(parse(text = core_code, keep.source = TRUE), envir = core)

core$selected_k <- c(27L, 28L, 38L, 44L)
evalq(
  build_partitions <- function(model, view) {
    fit <- stats::hclust(
      stats::dist(view$map$codes),
      method = "ward.D2"
    )
    partitions <- list()
    metadata <- list()
    for (k in selected_k) {
      solution_id <- paste0("inflect_ward_k", k)
      partitions[[solution_id]] <- validate_partition(
        stats::cutree(fit, k = k),
        k,
        view$map$nNodes,
        solution_id
      )
      metadata[[solution_id]] <- data.frame(
        solution_id = solution_id,
        family = "INFLECT Ward.D2 selected-k supplement",
        k = k,
        code_space = "model effective code matrix",
        seed = NA_integer_,
        stringsAsFactors = FALSE
      )
    }
    list(
      partitions = partitions,
      solution_metadata = do.call(rbind, metadata),
      inflect_codes = view$map$codes,
      historical_codes = NULL,
      consensus_ari = data.frame(
        solution_a = character(),
        solution_b = character(),
        ari = numeric(),
        stringsAsFactors = FALSE
      ),
      consensus_stable = NA,
      display_controls = character()
    )
  },
  envir = core
)

core$wrapper_file <- wrapper_file
core$core_file <- core_file
evalq(
  source_provenance <- function() {
    source_files <- c(
      wrapper_file,
      core_file,
      list.files(
        file.path(repo_root, "R"),
        pattern = "[.]R$",
        full.names = TRUE
      ),
      file.path(repo_root, "DESCRIPTION")
    )
    git_output <- function(...) {
      tryCatch(
        system2(
          "git",
          c("-C", repo_root, ...),
          stdout = TRUE,
          stderr = FALSE
        ),
        error = function(e) NA_character_
      )
    }
    list(
      git_commit = git_output("rev-parse", "HEAD")[[1L]],
      git_status = git_output("status", "--short"),
      source_md5 = unname(tools::md5sum(source_files)),
      source_files = normalizePath(source_files, mustWork = TRUE),
      selected_k = selected_k,
      implementation = paste0(
        "selected-k wrapper around the frozen primary modality workflow"
      )
    )
  },
  envir = core
)

switch(
  core$mode,
  stage = core$run_stage(),
  pilot = core$run_pilot(),
  base = core$run_base(),
  classify = core$run_classify(),
  refine = core$run_refine(),
  summarize = core$run_summarize(),
  status = core$status_report(),
  `self-test` = core$run_self_test(),
  stop(
    "Unknown mode `",
    core$mode,
    "`. Use stage, pilot, base, classify, refine, summarize, status, or self-test.",
    call. = FALSE
  )
)
