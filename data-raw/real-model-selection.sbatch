#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --output=data-raw/.real-model-selection-work/selection-%j.out
#SBATCH --error=data-raw/.real-model-selection-work/selection-%j.err

set -euo pipefail

repo_root=${SLURM_SUBMIT_DIR:?Submit this job from the fastINFLECT repository}
rscript_bin=${INFLECT_SELECTION_RSCRIPT:-/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export INFLECT_SELECTION_WORKDIR=${INFLECT_SELECTION_WORKDIR:-${repo_root}/data-raw/.real-model-selection-work}
export INFLECT_SELECTION_EVIDENCE_DIR=${INFLECT_SELECTION_EVIDENCE_DIR:-${repo_root}/inst/benchmarks/real-model-selection}
export INFLECT_MODALITY_WORKDIR=${INFLECT_MODALITY_WORKDIR:-${repo_root}/data-raw/.real-model-modality-work}

mkdir -p "${INFLECT_SELECTION_WORKDIR}"
cd "${repo_root}"
exec "${rscript_bin}" \
  data-raw/evaluate-real-model-selection.R \
  --mode=run
