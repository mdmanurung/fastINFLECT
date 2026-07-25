# Migrating to fastINFLECT 2.0

fastINFLECT 2.0 is a breaking release. It makes the tested cluster-count
schedule literal, retains complete transformed marker distributions by
default, separates dip and IQR evidence, and distinguishes fitted
estimates from partitions that actually exist.

The aggregate diagnostic remains useful for screening candidate values
of `k`. It does not establish that a cluster is truly unimodal.

## Supply a literal schedule

`set.i` is required. It must contain at least five unique, strictly
increasing integer values, and no value may exceed the number of SOM
nodes.

``` r
schedule <- seq.int(25L, 100L, by = 5L)
print(schedule)
result <- INFLECT(model, set.i = schedule)
```

The call evaluates exactly `25, 30, ..., 100`. The former shorthand
`set.i = c(25, 100)` now fails because its spacing is ambiguous.

For the pre-2.0 adaptive spacing, construct the schedule explicitly:

``` r
inflect_adaptive_set_i(n_nodes = 900L, max_k = 100L)
#>  [1]   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23
#> [20]  24  25  30  35  40  45  50  55  60  65  70  75  80  85  90  95 100
```

`max_k` is a hard upper bound. The helper performs overflow-safe
arithmetic and refuses schedules that are unsafe to materialise.

## Retain complete transformed values

`zeroes.in = TRUE` is the new default. It retains negative, zero, and
positive finite values. This matters for centred or otherwise
transformed marker data: dropping non-positive measurements changes the
tested distribution.

``` r
result <- INFLECT(
  model,
  set.i = schedule,
  zeroes.in = TRUE
)
```

`zeroes.in = FALSE` remains available as an explicit compatibility
setting. It excludes every non-positive value, warns when negative
inputs are present, and records exclusions:

``` r
compatibility_result <- INFLECT(
  model,
  set.i = schedule,
  zeroes.in = FALSE
)

compatibility_result$provenance$zero_handling
```

No zero-valued padding observations are manufactured.

## Read separated criteria

`uniform.test` selects the aggregate pass-rate curve but no longer hides
the component evidence:

| `uniform.test`  | Aggregate curve  | Evidence still retained |
| --------------- | ---------------- | ----------------------- |
| `"both"`        | dip and IQR pass | dip, IQR, combined      |
| `"unimodality"` | dip pass         | dip, IQR, combined      |
| `"spread"`      | IQR pass         | dip, IQR, combined      |

``` r
result <- INFLECT(
  model,
  set.i = schedule,
  uniform.test = "both"
)

result$scores
result$dip_pass[["50"]]
result$iqr_pass[["50"]]
result$combined_pass[["50"]]

details <- result$qc.details[["50"]]
details$dip_p_value
details$iqr
details$event_count
details$test_event_count
details$failure_reason
```

The canonical aggregate column is `qc_pass_rate`. The following aliases
are deprecated but retained during migration:

| Deprecated                           | Canonical        |
| ------------------------------------ | ---------------- |
| `U.set`, `collection.U`              | `scores`         |
| `Unimodality`                        | `qc_pass_rate`   |
| `Accuracy.sets`, `Accuracy.matrixes` | `criterion_pass` |

## Opt in to parallel scoring

Serial scoring is the default. Parallel scoring requires at least two
workers:

``` r
result <- INFLECT(
  model,
  set.i = schedule,
  multicore = TRUE,
  cores = 2L
)
```

Unix-like systems use forked workers over distinct SOM-node subtrees.
Unsupported platforms validate the same inputs and use a recorded serial
fallback. Every worker result is checked for schema, marker order, and
length before matrices are assembled.

Simulated dip p-values and optional event sampling use deterministic,
overflow-safe streams derived from the subtree and marker. They preserve
the caller’s RNG state and return the same result with one or two
workers.

## Bound event retention only as a sensitivity analysis

`max.events.per.node` samples each node once with `seed`. The same
retained rows feed every marker, both criteria, every k, and every
worker:

``` r
result <- INFLECT(
  model,
  set.i = schedule,
  max.events.per.node = 2000L,
  seed = 42L
)

result$provenance$scoring_mode
result$provenance$n_events_retained
result$provenance$event_sampling$per_node
```

The default `NULL` retains every event. A cap should only be adopted
after checking score and endpoint stability across relevant sample sizes
and seeds.

## Distinguish estimates from partitions

`result$selection` reports `partition_available`, `k_status`, and
`score_source`:

``` r
result$selection[c(
  "method",
  "k",
  "qc_pass_rate_at_k",
  "partition_available",
  "k_status",
  "score_source"
)]
```

A fitted-curve inflection can be an unscheduled integer. In that case it
is labelled `fitted_estimate_no_partition`; fastINFLECT does not
silently evaluate it. To obtain the partition and its criterion
evidence, include that literal k in a new schedule.

No value of k should be promoted solely because it sits at the
inflection of an aggregate pass-rate curve.

## Audit provenance before interpretation

The QC boundary records the criterion, thresholds, marker set, zero
handling, sampling order, seeds, worker backend, model layer weights,
retained events, and per-criterion summaries:

``` r
str(result$provenance, max.level = 2)
```

Use the following language for downstream validation:

  - `detected_multimodality` when independent tests reproducibly reject;
  - `no_detected_multimodality` when valid tests do not reject;
  - `ambiguous` for test disagreement, sensitivity, or borderline
    evidence;
  - `unresolved_discrete` for insufficient unique values, excessive
    ties, saturation, or test failure.

“No detected multimodality” is deliberately weaker than “unimodal.” The
dip test can reject a unimodal null; failure to reject does not prove
that null.
