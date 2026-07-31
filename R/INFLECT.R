#' @title Run a fastINFLECT metaclustering scan
#'
#' @description Evaluates every cluster count supplied in `set.i`, scores each
#' cluster-marker pair with separate dip and IQR criteria, and reports candidate
#' values of `k`. The result retains every tested partition, the component QC
#' evidence, and the settings needed to interpret the run. A candidate is a
#' screening result, not proof that a partition is biologically valid.
#'
#' @param FlowSOM.results A supported SOM object with completed SOM clustering. Supports \pkg{FlowSOM} objects and \pkg{kohonen} objects returned by \code{\link[kohonen]{som}} or \code{\link[kohonen]{xyf}}.
#' @param set.i Required vector of at least five unique, strictly increasing
#' integer cluster counts within the SOM-node range. Values are used literally.
#' See \code{\link{inflect_adaptive_set_i}} to construct an explicitly bounded
#' adaptive schedule.
#' @param multicore Use parallel QC scoring on platforms with fork support.
#'   Default `FALSE`.
#' @param cores Worker count when `multicore = TRUE`. Must be at least two. A
#'   `NULL` value uses the detected core count minus one.
#' @param zeroes.in If `TRUE` (default), retain negative, zero, and positive
#' finite transformed values. If `FALSE`, exclude all non-positive values,
#' record per-marker counts, and warn when negatives are present.
#' @param only.clustering.markers Evaluate only clustering markers. For a
#'   kohonen object, these are the markers in the first data layer.
#' @param acquired_markers Marker names to evaluate when
#'   `only.clustering.markers = FALSE`.
#' @param basedata Use the fitted `"Curve"` or observed `"Points"` to estimate
#'   the inflection.
#' @param ggtitle Optional title for the diagnostic plot.
#' @param uniform.test Criterion used for the aggregate pass-rate curve:
#'   `"both"` requires dip and IQR to pass, `"unimodality"` uses dip, and
#'   `"spread"` uses IQR. Component evidence is retained in every run.
#' @param th.pvalue Dip-test pass threshold. A pair passes the dip criterion
#'   when its p-value is at least this value. Default is \code{0.05}.
#' @param th.IQR Threshold for rejecting marker distribution based on inter-quartile range. Default is arc-sinh transformed value of \code{2}.
#' @param max.n.diptest Optional integer cap on the number of events used per
#'   cluster-marker dip test. The dip test's power depends on sample size, so a
#'   cap can be used for a sensitivity analysis. Subsampling is seeded and
#'   deterministic; the default \code{NULL} retains all testable events.
#' @param max.events.per.node Optional positive integer. If set, each SOM node
#' is sampled once to at most this many events. That shared event sample is used
#' for both dip and IQR scoring across all metaclusterings and workers. Default
#' \code{NULL} uses all events.
#' @param seed Non-negative integer base seed for optional event sampling.
#' Default \code{1L}.
#' @param target Target QC pass rate (fraction in \code{(0,1]} or percentage in
#' \code{(1,100]}) used to report the smallest tested k reaching the threshold.
#' Default \code{0.95}. See \code{\link{inflect_threshold_k}}.
#' @param verbose Print progress messages. Default `FALSE`.
#' @param ... Arguments to pass to \code{\link[diptest]{dip.test}} through \code{\link{FlowSOMQC}}.
#'
#' @return An S3 `inflect.results` object. Use `print()` for candidate values,
#'   `plot()` for the pass-rate curve, `as.data.frame()` for tested scores,
#'   `result$selection` for candidate metadata, and the criterion matrices plus
#'   `qc.details` for marker-level interpretation.
#'
#' @seealso \code{\link{iteration.metacluster}}, \code{\link{iteration.QC}}, \code{\link{FlowSOMQC}}, \code{\link{QC.to.curve}}, \code{\link{leastError}}, \code{\link{Lfunction}}, \code{\link{marker.performance}}
#'
#' @examples
#'
#' # Load the bundled, downsampled Levine32 FlowSOM object.
#' flowsom <- system.file(
#'   "extdata",
#'   "Levine32sample.Rdata",
#'   package = "fastINFLECT"
#' )
#' load(flowsom)
#'
#' result <- INFLECT(
#'   FlowSOM.results = dataset,
#'   set.i = 5:12,
#'   uniform.test = "both"
#' )
#'
#' result
#' plot(result)
#' as.data.frame(result)
#' result$selection
#'
#' @export
INFLECT <-
  function(FlowSOM.results,
           set.i,
           multicore = FALSE,
           cores = NULL,
           zeroes.in = TRUE,
           only.clustering.markers = TRUE,
           acquired_markers = NULL,
           basedata = "Curve",
           ggtitle = NULL,
           uniform.test = c("both", "spread", "unimodality"),
           th.pvalue = 0.05,
           th.IQR = 2,
           verbose = FALSE,
           max.n.diptest = NULL,
           max.events.per.node = NULL,
           seed = 1L,
           target = 0.95,
           ...) {
    if (missing(set.i)) {
      stop(
        "`set.i` is required; supply at least five literal cluster counts.",
        call. = FALSE
      )
    }
    start_time <- proc.time()[["elapsed"]]
    diptest_args <- list(...)
    uniform.test <- match.arg(uniform.test)
    invisible(.inflect_target_percent(target))
    cores <- .inflect_resolve_cores(multicore, cores)
    FlowSOM.results <- as_inflect_som(FlowSOM.results)

    set.i <- normalize_set_i(set.i, FlowSOM.results)
    markers <- resolve_inflect_markers(
      FlowSOM.results = FlowSOM.results,
      only.clustering.markers = only.clustering.markers,
      acquired_markers = acquired_markers
    )

    metaclustering.list <-
      iteration.metacluster(
        FlowSOM.results = FlowSOM.results,
        set.i = set.i,
        multicore = multicore,
        cores = cores
      )

    qc <-
      iteration.QC(
        FlowSOM.results = FlowSOM.results,
        metaclustering.list = metaclustering.list,
        set.i = set.i,
        multicore = multicore,
        cores = cores,
        zeroes.in = zeroes.in,
        only.clustering.markers = only.clustering.markers,
        acquired_markers = acquired_markers,
        uniform.test = uniform.test,
        th.pvalue = th.pvalue,
        th.IQR = th.IQR,
        max.n.diptest = max.n.diptest,
        max.events.per.node = max.events.per.node,
        seed = seed,
        verbose = verbose,
        ...
      )


    diagnostic.graph <-
      QC.to.curve(collection.U = qc,
                  basedata = basedata,
                  ggtitle = ggtitle)
    diagnostic_scores <- diagnostic.graph$scores
    if (is.null(diagnostic_scores)) {
      diagnostic_scores <- .inflect_score_frame(diagnostic.graph$collection.U)
    }

    selection <- .inflect_selection(
      collection.U = diagnostic_scores,
      lfunction = diagnostic.graph$lfunction,
      fittedcurve = diagnostic.graph$fittedcurve,
      target = target
    )

    provenance <- build_inflect_provenance(
      FlowSOM.results = FlowSOM.results,
      set.i = set.i,
      uniform.test = uniform.test,
      th.pvalue = th.pvalue,
      th.IQR = th.IQR,
      zeroes.in = zeroes.in,
      basedata = basedata,
      only.clustering.markers = only.clustering.markers,
      acquired_markers = acquired_markers,
      markers = markers,
      elapsed_seconds = proc.time()[["elapsed"]] - start_time,
      diptest_args = diptest_args,
      max.n.diptest = max.n.diptest,
      max.events.per.node = max.events.per.node,
      seed = seed,
      target = target,
      qc_provenance = qc$provenance
    )

    return(new_inflect_results(
      scores = diagnostic_scores,
      collection.U = diagnostic.graph$collection.U,
      fittedcurve = diagnostic.graph$fittedcurve,
      lfunction = diagnostic.graph$lfunction,
      ggplot = diagnostic.graph$ggplot,
      metaclustering.list = metaclustering.list,
      accuracy.sets = if (is.null(qc$criterion_pass)) {
        qc$Accuracy.matrixes
      } else {
        qc$criterion_pass
      },
      qc.details = if (is.null(qc$qc.details)) list() else qc$qc.details,
      dip_pass = if (is.null(qc$dip_pass)) list() else qc$dip_pass,
      iqr_pass = if (is.null(qc$iqr_pass)) list() else qc$iqr_pass,
      combined_pass = if (is.null(qc$combined_pass)) list() else qc$combined_pass,
      criterion_pass = if (is.null(qc$criterion_pass)) {
        qc$Accuracy.matrixes
      } else {
        qc$criterion_pass
      },
      selection = selection,
      provenance = provenance
    ))
  }
