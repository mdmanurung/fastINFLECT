# Large-model per-node cap stability benchmark

> **Invalid as modality evidence.** This historical engineering run used
> `uniform.test = "spread"` (IQR only), excluded all non-positive transformed
> values, and did not materialise its fitted k = 62 estimate. It can support the
> recorded runtime and cap-sensitivity observations only. It cannot support a
> claim of unimodality or “no detected multimodality.” The replacement audit is
> `data-raw/validate-real-model-modality.R`.

Date: 2026-07-24
Base commit: `2390ae2974131beccd47950a95a165a8f1b0cf1a` plus the fastINFLECT
2.0 working-tree implementation described below.

## Contract

- Model: 39,050,953 events, 900 SOM nodes, 27 data columns; source RDS size
  7.9 GiB.
- Schedule: `seq.int(25L, 100L, by = 5L)`; no k above 100 was evaluated.
- Scoring: `uniform.test = "spread"`, `zeroes.in = FALSE`.
- Execution: `multicore = TRUE`, `cores = 2L`.
- Candidate `max.events.per.node`: 1,000, 2,000, and 5,000.
- Candidate seeds: 1, 42, and 2026.
- Acceptance: at most 1 percentage point maximum score deviation from the
  uncapped run; every recommended k within 5 of uncapped; stable across seeds.

The model loaded in 90.9 seconds. The complete benchmark took 6 minutes
43.57 seconds and had a maximum resident set size of 12,464,840 KiB
(approximately 11.89 GiB) according to `/usr/bin/time -v`. The uncapped
two-worker scoring pass took 160.753 seconds and completed without an OOM or
worker-result error.

## Results

| Cap | Seed | Retained events | Wall time (s) | Max score deviation (pp) | Inflection k | Kneedle k | Max recommendation delta | Pass |
|---:|---:|---:|---:|---:|---:|---:|---:|:---:|
| uncapped | 1 | 39,050,953 | 160.753 | 0.000 | 62 | 55 | 0 | yes |
| 1,000 | 1 | 897,839 | 9.238 | 2.074 | 32 | 55 | 30 | no |
| 1,000 | 42 | 897,839 | 7.965 | 2.074 | 33 | 55 | 29 | no |
| 1,000 | 2026 | 897,839 | 8.240 | 1.778 | 29 | 50 | 33 | no |
| 2,000 | 1 | 1,793,851 | 12.225 | 2.222 | 29 | 50 | 33 | no |
| 2,000 | 42 | 1,793,851 | 11.972 | 1.926 | 33 | 55 | 29 | no |
| 2,000 | 2026 | 1,793,851 | 12.012 | 2.074 | 29 | 50 | 33 | no |
| 5,000 | 1 | 4,479,461 | 23.969 | 2.222 | 29 | 50 | 33 | no |
| 5,000 | 42 | 4,479,461 | 23.928 | 2.222 | 30 | 50 | 32 | no |
| 5,000 | 2026 | 4,479,461 | 24.048 | 2.074 | 28 | 50 | 34 | no |

The 95% threshold recommendation was `NA` in every run. The uncapped
inflection k of 62 was fitted rather than directly tested; kneedle 55 was
directly tested.

All candidates were internally stable across seeds: their maximum score ranges
were 0.741, 0.494, and 0.212 percentage points for caps 1,000, 2,000, and
5,000, respectively, and recommendation ranges were at most one five-k step.
Nevertheless, every cap failed comparison with the uncapped source of truth:
score deviations exceeded 1 percentage point and fitted inflection k shifted
by 29 to 34.

Decision: leave `max.events.per.node = NULL` for this notebook and model.

## Exported-call notebook acceptance

The final notebook call was also run end to end through exported `INFLECT()`:

- `set.i`: 25, 30, ..., 100
- `multicore = TRUE`, `cores = 2L`
- `uniform.test = "spread"`
- `max.events.per.node = NULL`
- `seed = 42L`
- `INFLECT()` wall time: 168.958 seconds
- Whole process including model loading: 4 minutes 25.15 seconds
- Whole-process maximum resident set size: 17,387,288 KiB
  (approximately 16.58 GiB)
- Original/retained events: 39,050,953 / 39,050,953
- Recommendation: fitted inflection k 62 at 73.63536%; directly tested
  kneedle 55 at 74.07407%; 95% threshold not reached
- Warnings signaled inside `INFLECT()`: none

The run completed with no k above 100, interpolation-warning flood, OOM,
worker failure, or downstream `dimnames` error.

## Artifacts and source identity

- Machine-readable result:
  `inst/benchmarks/large-model-cap-stability-2026-07-24.rds`
- Exported-call acceptance:
  `inst/benchmarks/large-model-notebook-acceptance-2026-07-24.rds`
- Reproduction script: `data-raw/benchmark-large-model-cap.R`
- Source SHA-256:
  - `R/INFLECT.R`: `5a5ddab873ae4639469fab07b31b25998bb3156d34f0c61b1f6ed38ff45b1907`
  - `R/iteration-QC.R`: `5c813979dad77623e94c14c638db96359356e847a7681e59ec9d2ddc82c412e1`
  - `R/inflect-qc-core.R`: `d81ff63e8ce804a6c3836176899ed86395a93ad7c4d8f45adfc363ab3183821c`
  - `R/inflect-provenance.R`: `e1ed1cc2b42ccab22cce3acf86eb7430235f566172a9822dd40d60df52757360`
  - `R/select-k.R`: `aa623f5f02a572f04d2665b0245ef9a0b844424bcd0cfcb68e6337897a0066e3`
  - `data-raw/benchmark-large-model-cap.R`:
    `66192fe7c6340eee621456766deb4afdbec77bf8644d74b48ecd89d4ec0d185a`
