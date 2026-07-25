# Historical large-model artifacts

The following 2026-07-24 artifacts are retained for provenance but are invalid
as evidence of cluster modality:

- `large-model-cap-stability-2026-07-24.rds`
- `large-model-notebook-acceptance-2026-07-24.rds`
- `large-model-cap-stability-2026-07-24.md`

They used the IQR-only `uniform.test = "spread"` curve,
`zeroes.in = FALSE`, and a sparse schedule that did not materialise the fitted
k = 62 estimate. The results may be cited for runtime, memory, and
cap-sensitivity observations only.

## Valid replacement

The replacement workflows completed on 2026-07-25:

- `real-model-selection/selection-evidence.rds` contains the literal 25:100
  INFLECT sweep, repaired event-weighted metrics, three-seed consensus
  comparison, and directly materialised candidate partitions.
- `selected-k-modality/modality-evidence.rds` directly audits k = 27, 28, 38,
  and 44.
- `real-model-modality/modality-evidence.rds` audits k = 55, 60, 62, and 65
  plus historical 80/20 Ward and seeded FlowSOM-consensus comparators.
- `real-model-selection/bmv-k-selection-synthesis.rds` joins those independent
  artifacts and records the operational recommendation and sensitivity range.

Both modality bundles passed their completeness gates: every requested
partition was materialised, every base and refinement combination was
accounted for, and no interpolation warning was detected. The valid evidence
supports an operational k recommendation but not a claim of global
unimodality.
