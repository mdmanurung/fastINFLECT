# Real-model modality audit

`validate-real-model-modality.R` is the resumable scientific audit for the
39,050,953-event EXV BMV model. Its scope is univariate marker marginals only.

## Analysis dependency

Ordinary fastINFLECT use does not require `multimode`. Install it into the
ignored task-local library for this audit:

```bash
mkdir -p data-raw/.real-model-lib
R_MAKEVARS_USER=data-raw/Makevars.audit \
  /exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript \
  -e "install.packages('multimode', lib='data-raw/.real-model-lib', repos='https://cloud.r-project.org')"
```

The workflow uses `multimode::modetest(mod0 = 1, method = "ACR", B = 1999)`.
See the [multimode reference
page](https://search.r-project.org/CRAN/refmans/multimode/html/modetest.html)
and the [ACR methodology](https://arxiv.org/abs/1609.05188).

## Stages

The stages are dependency ordered:

1. `stage` loads the model; materialises INFLECT Ward.D2 k = 55, 60, 62, 65,
   historical intended 80/20 Ward.D2 k = 50 and 62, and FlowSOM consensus
   k = 62 for seeds 1, 42, and 2026; then writes nested 5,000-event sample
   checkpoints.
2. `pilot` runs exactly 100 ACR tests at B = 1999 and projects base wall time.
3. `base` runs the 2,000-event tests in zero-based array shards.
4. `classify` applies BH correction within each 27-marker
   solution-cluster-seed-size family and writes the refinement manifest.
5. `refine` repeats required clusters at 1,000 and 5,000 events.
6. `summarize` checks completeness, integrates sample-size sensitivity, and
   writes the evidence bundle and plots.

Run the complete dependency graph with:

```bash
INFLECT_ARRAY_COUNT=64 INFLECT_ARRAY_LIMIT=16 \
  data-raw/submit-real-model-modality.sh
```

The array count controls deterministic task sharding; it does not change
samples or test seeds. Existing per-cluster checkpoints are reused.

On schedulers whose maximum array index is 63, the optional
`real-model-modality-reshard.sbatch` wrapper maps two physical 0:63 arrays onto
odd logical shards 1:127 and 129:255. This is only a checkpoint-resume
acceleration for the slower 5,000-event refinement rows; it does not alter
samples, tests, seeds, or output paths.

## Cluster-number comparison

`evaluate-real-model-selection.R` complements the modality audit with the
criteria needed to compare candidate values of k. It requires finite marker
values, retains the complete full-event distributions, runs the literal
`25:100` schedule, and keeps separate dip, IQR-spread, and combined pass rates.
It also compares:

- silhouette, Calinski-Harabasz, Davies-Bouldin, and explained-dispersion
  knees on both the fitted model code space and the historical intended 80/20
  code space;
- FlowSOM/ConsensusClusterPlus CDF area and seed-to-seed adjusted Rand index
  across k;
- deterministic 5,000-event sampled IQR evidence for every Ward and consensus
  control materialised by the modality stage.

Run it after the modality `stage` job has completed:

```bash
sbatch --dependency=afterok:<stage-job-id> \
  --time=12:00:00 --mem=96G \
  data-raw/real-model-selection.sbatch
```

The output is checkpointed under ignored
`data-raw/.real-model-selection-work/` and written to
`inst/benchmarks/real-model-selection/`. Consensus delta area is reported as a
heuristic plateau, not as an automatic proof that one k is biologically
correct.

Each 100-replicate consensus seed is slow but checkpoint-separable. The BMV
seed-2026 helper used about 1.84 GiB peak RSS; the near-96-GiB peak belongs to
the earlier full-event INFLECT/fork stage, not to consensus. An untouched seed
can be run independently, with identical numerical parameters, using:

```bash
sbatch --time=12:00:00 --mem=8G \
  data-raw/real-model-consensus-seed.sbatch <seed>
```

The helper writes the same compact seed checkpoint and additionally records
its script, input-checkpoint, parameter, dependency, and session provenance.
Do not launch it for a seed already being computed by the serial job.

Event-weighted dispersion excludes zero-event SOM nodes from centroid
calculations while retaining them in nominal-k and empty-cluster diagnostics.
This avoids undefined zero-weight centroids without hiding empty clusters.
When source or metric logic changes, use a new
`INFLECT_SELECTION_WORKDIR`; do not reuse a stale run contract.

If that dense comparison selects values of k that were not included in the
primary modality candidates, run the selected-k supplement:

```bash
INFLECT_ARRAY_COUNT=64 INFLECT_ARRAY_LIMIT=64 \
  data-raw/submit-selected-k-modality.sh
```

The supplement uses the same frozen sampling, dip, ACR, BH, classification,
and refinement implementation. Its current declared scope is k = 27, 28, 38,
and 44; k = 55 is already present in the primary audit. Outputs are written to
`inst/benchmarks/selected-k-modality/`.

Once selection and both modality bundles are complete, run:

```bash
/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript \
  data-raw/summarize-bmv-k-selection.R
```

This fail-closed synthesis checks model identity, completeness, finite repaired
metrics, and direct selected-k coverage before writing the operational
recommendation, exact comparison tables, and provenance bundle.

The dense LL.4 fits can emit `NaNs produced` while the optimiser explores a
temporary nonpositive ED50. Run
`data-raw/diagnose-bmv-curve-warning.R` to verify the warning call, fitted-value
finiteness, and exact reproduction of all stored inflections. The diagnostic
retains the warning; it does not suppress it.

After both final evidence bundles exist, execute the BMV notebook with:

```bash
sbatch --dependency=afterok:<modality-summary-job-id>:<selected-k-summary-job-id>:<selection-job-id> \
  --time=12:00:00 --mem=96G \
  data-raw/real-model-notebook.sbatch
```

The notebook job writes to an ignored temporary path, verifies that every code
cell executed, and only then replaces `data-raw/test_inflect.ipynb`.

Inspect progress without loading the model:

```bash
/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript \
  data-raw/validate-real-model-modality.R --mode=status
```

Run the lightweight workflow-contract test with `--mode=self-test`.

## Evidence contract

Every requested partition must have exactly k labels and cover all 900 nodes.
Every cluster-marker-seed combination must be present or carry an explicit
failure reason. Tied or saturated measurements are never jittered.

The base status rules are:

- `detected_multimodality`: both dip and ACR reject at q < 0.05 in at least
  two seeds;
- `no_detected_multimodality`: neither rejects in any seed, all runs are
  valid, and no q lies in [0.05, 0.10);
- `ambiguous`: discordant tests, seed instability, borderline q, or
  non-replicated rejection;
- `unresolved_discrete`: insufficient unique values, excessive ties, boundary
  saturation, or test failure.

Detected and ambiguous pairs plus a deterministic, marker- and
cluster-size-stratified 5% audit of apparent passes are repeated at 1,000 and
5,000 events. Sample-size disagreement is integrated conservatively as
ambiguous; an unresolved required refinement remains unresolved.

BH correction controls the 27-marker family within a cluster, method, seed,
and sample size. It is not global correction across clusters.

## Outputs

Checkpoints live under the ignored
`data-raw/.real-model-modality-work/`. Final evidence is written to
`inst/benchmarks/real-model-modality/`, including:

- the complete RDS bundle;
- base and refinement test CSVs;
- base and final classifications;
- refinement sensitivity;
- partition labels and node-level adjusted Rand indices;
- weighted and unweighted status summaries;
- solution statements;
- status heatmaps and distribution audits;
- recorded source/model identity, dependency versions, peak RSS, and wall
  time.

A global “no detected multimodality” statement is allowed only when every
assessed marker in every cluster has concordant valid evidence with no
detected, ambiguous, or unresolved entry. It must never be restated as proof
that the solution is truly unimodal.
