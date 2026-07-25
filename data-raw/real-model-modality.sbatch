#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --output=data-raw/.real-model-modality-work/logs/%x-%A-%a.out
#SBATCH --error=data-raw/.real-model-modality-work/logs/%x-%A-%a.err

set -euo pipefail

repo_root=${SLURM_SUBMIT_DIR:?Submit this job from the fastINFLECT repository}
mode=${INFLECT_MODALITY_MODE:?Set INFLECT_MODALITY_MODE}
rscript_bin=${INFLECT_MODALITY_RSCRIPT:-/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export INFLECT_MODALITY_WORKDIR=${INFLECT_MODALITY_WORKDIR:-${repo_root}/data-raw/.real-model-modality-work}
export INFLECT_MODALITY_EVIDENCE_DIR=${INFLECT_MODALITY_EVIDENCE_DIR:-${repo_root}/inst/benchmarks/real-model-modality}
export INFLECT_MODALITY_R_LIB=${INFLECT_MODALITY_R_LIB:-${repo_root}/data-raw/.real-model-lib}

cd "${repo_root}"
exec "${rscript_bin}" \
  data-raw/validate-real-model-modality.R \
  "--mode=${mode}"
