# fastINFLECT

fastINFLECT evaluates a literal set of metacluster counts for a FlowSOM
or kohonen self-organising map. It returns the tested partitions,
marker-level QC evidence, and candidate values of `k`. Use the
candidates to narrow the search, then inspect the retained evidence
before choosing a partition.

## Install fastINFLECT

``` r

# install.packages("pak")
pak::pak("mdmanurung/fastINFLECT")
```

Install `kohonen` only when the input is a kohonen SOM:

``` r

install.packages("kohonen")
```

## Run an analysis

The package includes a downsampled Levine32 FlowSOM object.

``` r

library(fastINFLECT)

load(system.file(
  "extdata",
  "Levine32sample.Rdata",
  package = "fastINFLECT"
))

result <- INFLECT(
  FlowSOM.results = dataset,
  set.i = 5:12,
  uniform.test = "both"
)
```

`set.i` is the exact schedule to evaluate. Supply at least five unique,
increasing integers within the number of SOM nodes. The default
`zeroes.in = TRUE` retains every finite transformed value, including
zero and negative values.

## Read the result

``` r

result
plot(result)
as.data.frame(result)
result$selection
```

The main fields answer different questions:

| Field | Meaning |
|----|----|
| `scores` | QC pass rate for each tested `k` |
| `selection` | Inflection, Kneedle, and target-threshold candidates |
| `metaclustering.list` | Node labels for every tested partition |
| `criterion_pass` | Decision matrix selected by `uniform.test` and used for `qc_pass_rate` |
| `dip_pass`, `iqr_pass`, `combined_pass` | Separate cluster-by-marker decisions |
| `qc.details` | P-values, IQRs, event counts, exclusions, and failure reasons |
| `provenance` | Schedule, markers, thresholds, sampling, and run settings |

Candidate rows are diagnostic summaries, not a ranking. The example
below uses the directly tested inflection candidate only to demonstrate
inspection; it is not a default recommendation. Only a row with
`partition_available = TRUE` has a matching partition. If a finite
candidate was not directly tested, inspect `k_status`, add its `k` to
`set.i`, and rerun before inspection. Candidate methods that happen to
return the same `k` are not independent confirmation.

## Inspect a candidate partition

``` r

candidate_method <- "inflection"
candidate <- result$selection[
  result$selection$method == candidate_method,
  ,
  drop = FALSE
]

if (nrow(candidate) != 1L) {
  stop("Expected exactly one inflection candidate.", call. = FALSE)
}
if (!isTRUE(candidate$partition_available[[1L]])) {
  stop(
    "Inspect `k_status`; if the candidate k is finite, add it to `set.i` and rerun.",
    call. = FALSE
  )
}

candidate_k <- candidate$k[[1L]]
key <- as.character(candidate_k)
if (!key %in% names(result$metaclustering.list)) {
  stop("The directly tested candidate partition is missing.", call. = FALSE)
}

candidate[c("method", "k", "qc_pass_rate_at_k", "partition_available", "k_status")]

partition <- result$metaclustering.list[[key]]
partition_map <- data.frame(
  som_node = seq_along(partition),
  metacluster = unname(partition)
)
head(partition_map)
table(partition_map$metacluster)

criterion_summary <- result$provenance$criterion_summary
selected_criterion_summary <- criterion_summary[
  criterion_summary$k == candidate_k & criterion_summary$selected,
  c("k", "criterion", "passed", "failed", "unresolved", "total", "qc_pass_rate"),
  drop = FALSE
]
if (nrow(selected_criterion_summary) != 1L) {
  stop("Expected exactly one selected criterion summary.", call. = FALSE)
}
selected_criterion_summary

result$criterion_pass[[key]]
result$dip_pass[[key]]
result$iqr_pass[[key]]
result$combined_pass[[key]]
result$qc.details[[key]]$failure_reason

marker_qc <- marker.performance(result)
marker_qc$plot
```

This is an inspected candidate partition, not a scientifically selected
partition. Positions are SOM nodes and values are nominal metacluster
labels. The table counts SOM nodes, not events or cells. Labels are not
ordered scores and must not be compared numerically across different
values of `k`.

`criterion_pass[[key]]` is the matrix used to calculate `qc_pass_rate`.
It equals `dip_pass[[key]]` for `uniform.test = "unimodality"`,
`iqr_pass[[key]]` for `"spread"`, and `combined_pass[[key]]` for
`"both"`. `TRUE` passed the recorded threshold, `FALSE` did not, and
`NA` means that the selected decision was unavailable. Under `"both"`, a
definitive component failure can determine `FALSE` even if the other
component is unresolved, so inspect the component matrices and
`failure_reason` as well.

`qc_pass_rate` is `100 * passed / total` over the metacluster-by-marker
matrix selected by `uniform.test`; unresolved (`NA`) decisions remain in
`total`. Because `k`, cluster event counts, and test behavior change
across candidates, compare `passed`, `failed`, `unresolved`, and
`total`. A higher rate means only that a larger fraction passed. It does
not imply fewer absolute failures or that a partition is preferable or
biologically valid.

Use these interpretations:

- Dip-test non-rejection means no detected multimodality at the recorded
  threshold, not proof of unimodality.
- IQR is a spread screen, not a unimodality test.
- Per-pair thresholds are not multiplicity-controlled inference.
- Increasing `k` changes group sizes and test behavior, so a rising
  curve is diagnostic rather than evidence that larger `k` is
  biologically better.
- A candidate `k` remains diagnostic until its marker-level evidence,
  stability, and biological usefulness are checked.

The complete workflow is in
[`vignette("fastINFLECT")`](https://mdmanurung.github.io/fastINFLECT/articles/fastINFLECT.md).
