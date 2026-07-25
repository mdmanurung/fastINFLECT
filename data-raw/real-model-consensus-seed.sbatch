#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --output=data-raw/.real-model-selection-work/consensus-helper-%j.out
#SBATCH --error=data-raw/.real-model-selection-work/consensus-helper-%j.err

set -euo pipefail

if [[ $# -ne 1 ]] || ! [[ $1 =~ ^[1-9][0-9]*$ ]]; then
  echo "Usage: sbatch real-model-consensus-seed.sbatch <positive-seed>" >&2
  exit 2
fi

repo_root=${SLURM_SUBMIT_DIR:?Submit this job from the fastINFLECT repository}
rscript_bin=${INFLECT_SELECTION_RSCRIPT:-/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export INFLECT_SELECTION_WORKDIR=${INFLECT_SELECTION_WORKDIR:-${repo_root}/data-raw/.real-model-selection-work}

mkdir -p "${INFLECT_SELECTION_WORKDIR}"
cd "${repo_root}"
exec "${rscript_bin}" \
  data-raw/run-real-model-consensus-seed.R \
  "--seed=$1"
