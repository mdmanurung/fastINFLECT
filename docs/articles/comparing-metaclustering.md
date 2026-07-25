# Comparing metaclustering strategies with fastINFLECT

FlowSOM reduces events to SOM nodes; metaclustering assigns those nodes
to a smaller number of clusters. Different partition engines and
code-space weights can produce different node assignments even when they
use the same value of `k`. A fair comparison therefore separates three
questions:

1.  Which code matrix and weights were used?
2.  Which partition was actually materialised?
3.  What marker-level evidence remains inside that partition?

An aggregate fastINFLECT pass-rate curve answers none of these questions
by itself.

## Materialise each partition

The examples below are deliberately not executed while building the
vignette. They show the required contract for a user-supplied FlowSOM
object named `model`.

### Ward.D2 hierarchy

fastINFLECT uses a Ward.D2 hierarchy over the adapter’s effective code
matrix. The literal schedule is retained as named node-label vectors:

``` r
schedule <- c(25L, 30L, 35L, 40L, 45L, 50L, 55L, 60L, 62L, 65L)

ward_partitions <- iteration.metacluster(
  FlowSOM.results = model,
  set.i = schedule,
  multicore = FALSE
)

stopifnot(
  identical(names(ward_partitions), as.character(schedule)),
  all(vapply(
    seq_along(schedule),
    function(i) length(unique(ward_partitions[[i]])) == schedule[[i]],
    logical(1)
  ))
)
```

For a kohonen multi-layer SOM, the adapter multiplies each code layer by
the square root of its user weight and distance weight before
concatenating the layers. The recorded values are available in:

``` r
result$provenance$model
```

These weights are part of the partition definition. Equal fitted model
weights and an historical 80/20 construction are different analyses, not
interchangeable implementations.

### FlowSOM consensus control

FlowSOM consensus metaclustering is a separate partition engine:

``` r
consensus_seeds <- c(1L, 42L, 2026L)

consensus_partitions <- setNames(
  lapply(consensus_seeds, function(seed) {
    FlowSOM::metaClustering_consensus(
      data = effective_codes,
      k = 62L,
      seed = seed
    )
  }),
  paste0("seed_", consensus_seeds)
)
```

Compare all three node-label vectors. If every pairwise adjusted Rand
index is at least 0.95, seed 42 can be used as a displayed control while
retaining all three partitions in provenance. Otherwise, the controls
are seed-unstable and all three should be shown.

## Score a common marker panel

`iteration.QC()` accepts a named list containing exactly one partition
for every requested k. It returns separate dip, IQR, and combined
evidence:

``` r
qc <- iteration.QC(
  FlowSOM.results = model,
  metaclustering.list = ward_partitions,
  set.i = schedule,
  zeroes.in = TRUE,
  uniform.test = "both",
  multicore = TRUE,
  cores = 2L,
  seed = 42L
)

qc$scores
qc$dip_pass[["62"]]
qc$iqr_pass[["62"]]
qc$combined_pass[["62"]]
qc$qc.details[["62"]]$failure_reason
```

`uniform.test = "spread"` selects only the IQR pass rate for the
aggregate curve. It does not perform a dip-based modality assessment.
Conversely, `uniform.test = "unimodality"` selects the dip result but
still cannot prove a unimodal distribution.

The default `zeroes.in = TRUE` retains the complete finite transformed
distribution. Setting it to `FALSE` excludes every non-positive value
and can substantially change a centred marker distribution.

## Compare solutions without hiding cluster size

For each fully materialised solution, report at least:

  - unweighted fractions of cluster-marker statuses;
  - event-weighted fractions of the same statuses;
  - the fraction of clusters with no detected or unresolved marker;
  - pairwise node-level adjusted Rand index;
  - a cluster-by-marker status heatmap.

Event-weighted and unweighted summaries answer different questions. A
small problematic cluster counts equally in the unweighted summary but
contributes few events to the weighted summary. Both are needed.

## Validate modality independently

For consequential claims, fastINFLECT’s aggregate screen should be
followed by independent tests on deterministic samples from complete
transformed distributions. The repository workflow uses:

  - `diptest::dip.test()`;
  - `multimode::modetest(mod0 = 1, method = "ACR", B = 1999)`;
  - seeds 1, 42, and 2026;
  - Benjamini-Hochberg correction across the marker family within each
    cluster, method, seed, and sample size;
  - repeats at 1,000 and 5,000 events for detected or ambiguous pairs
    and for a deterministic 5% audit of apparent passes.

The correction controls the marker family within a cluster. It is not
global control across all clusters.

Tied or saturated cytometry measurements are not silently jittered.
Pairs with insufficient unique values, excessive ties, boundary
saturation, or test failure are labelled `unresolved_discrete`.

## Interpret conservatively

Use four outcome labels:

  - `detected_multimodality`: both independent tests reject after
    correction in at least two seeds and required sample-size checks
    agree;
  - `no_detected_multimodality`: neither test rejects in any valid seed
    and required checks remain stable;
  - `ambiguous`: test disagreement, seed instability, borderline
    evidence, or sample-size sensitivity;
  - `unresolved_discrete`: the measurement grid or a test failure
    prevents a valid assessment.

Do not call a cluster “truly unimodal.” Failure to reject a unimodal
null is not proof of that null, and univariate marker marginals do not
establish multivariate cluster homogeneity or biological validity.

The pre-2.0 `vignette_cache.rds` is retained only as historical package
performance material. Because it used `zeroes.in = FALSE` and an
aggregate pass-rate interpretation, it is not scientific evidence for
cluster modality.
