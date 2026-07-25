#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --output=data-raw/.real-model-selection-work-v2/selection-%j.out
#SBATCH --error=data-raw/.real-model-selection-work-v2/selection-%j.err

set -euo pipefail

repo_root=${SLURM_SUBMIT_DIR:?Submit this job from the fastINFLECT repository}
rscript_bin=${INFLECT_SELECTION_RSCRIPT:-/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript}
upstream_dir=${INFLECT_SELECTION_UPSTREAM_WORKDIR:-${repo_root}/data-raw/.real-model-selection-work}
target_dir=${INFLECT_SELECTION_WORKDIR:-${repo_root}/data-raw/.real-model-selection-work-v2}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export INFLECT_SELECTION_WORKDIR=${target_dir}
export INFLECT_SELECTION_EVIDENCE_DIR=${INFLECT_SELECTION_EVIDENCE_DIR:-${repo_root}/inst/benchmarks/real-model-selection}
export INFLECT_MODALITY_WORKDIR=${INFLECT_MODALITY_WORKDIR:-${repo_root}/data-raw/.real-model-modality-work}

mkdir -p "${target_dir}/consensus"
if [[ ! -s "${target_dir}/full-inflect-selection.rds" ]]; then
  cp \
    "${upstream_dir}/full-inflect-selection.rds" \
    "${target_dir}/full-inflect-selection.rds"
fi
for seed in 0001 0042 2026; do
  source_checkpoint="${upstream_dir}/consensus/consensus-seed-${seed}.rds"
  target_checkpoint="${target_dir}/consensus/consensus-seed-${seed}.rds"
  if [[ ! -s "${target_checkpoint}" ]]; then
    cp "${source_checkpoint}" "${target_checkpoint}"
  fi
  test -s "${target_checkpoint}"
done

cd "${repo_root}"
exec "${rscript_bin}" \
  data-raw/evaluate-real-model-selection.R \
  --mode=run
