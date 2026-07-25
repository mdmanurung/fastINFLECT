#!/usr/bin/env bash
#SBATCH --cpus-per-task=1
#SBATCH --output=data-raw/.real-model-modality-work/logs/%x-%A-%a.out
#SBATCH --error=data-raw/.real-model-modality-work/logs/%x-%A-%a.err

set -euo pipefail

repo_root=${SLURM_SUBMIT_DIR:?Submit this job from the fastINFLECT repository}
raw_shard=${SLURM_ARRAY_TASK_ID:?Run this wrapper as an array job}
logical_base=${INFLECT_RESHARD_LOGICAL_BASE:?Set logical shard base to 1 or 129}

if ! [[ ${raw_shard} =~ ^[0-9]+$ ]] ||
   ! [[ ${logical_base} =~ ^(1|129)$ ]]; then
  echo "Invalid reshard index or logical base." >&2
  exit 2
fi

logical_shard=$((2 * raw_shard + logical_base))
if (( logical_shard < 1 || logical_shard > 255 || logical_shard % 2 == 0 )); then
  echo "Computed logical shard is outside the odd 1:255 set." >&2
  exit 2
fi

export SLURM_ARRAY_TASK_ID=${logical_shard}
export INFLECT_ARRAY_COUNT=256
export INFLECT_MODALITY_MODE=refine

exec bash "${repo_root}/data-raw/real-model-modality.sbatch"
