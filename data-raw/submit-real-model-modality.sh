#!/usr/bin/env bash

set -euo pipefail

repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
array_count=${INFLECT_ARRAY_COUNT:-64}
array_limit=${INFLECT_ARRAY_LIMIT:-16}
job_script="${repo_root}/data-raw/real-model-modality.sbatch"
work_dir=${INFLECT_MODALITY_WORKDIR:-${repo_root}/data-raw/.real-model-modality-work}

if ! [[ ${array_count} =~ ^[1-9][0-9]*$ ]]; then
  echo "INFLECT_ARRAY_COUNT must be a positive integer." >&2
  exit 2
fi
if ! [[ ${array_limit} =~ ^[1-9][0-9]*$ ]]; then
  echo "INFLECT_ARRAY_LIMIT must be a positive integer." >&2
  exit 2
fi

mkdir -p "${work_dir}/logs"

job_id() {
  printf '%s' "${1%%;*}"
}

stage_raw=$(sbatch \
  --parsable \
  --chdir="${repo_root}" \
  --job-name=inflect-modality-stage \
  --time=12:00:00 \
  --mem=96G \
  --export=ALL,INFLECT_MODALITY_MODE=stage \
  "${job_script}")
stage_id=$(job_id "${stage_raw}")

pilot_raw=$(sbatch \
  --parsable \
  --chdir="${repo_root}" \
  --job-name=inflect-modality-pilot \
  --dependency="afterok:${stage_id}" \
  --time=06:00:00 \
  --mem=8G \
  --export=ALL,INFLECT_MODALITY_MODE=pilot \
  "${job_script}")
pilot_id=$(job_id "${pilot_raw}")

base_raw=$(sbatch \
  --parsable \
  --chdir="${repo_root}" \
  --job-name=inflect-modality-base \
  --dependency="afterok:${stage_id}" \
  --array="0-$((array_count - 1))%${array_limit}" \
  --time=12:00:00 \
  --mem=8G \
  --export="ALL,INFLECT_MODALITY_MODE=base,INFLECT_ARRAY_COUNT=${array_count}" \
  "${job_script}")
base_id=$(job_id "${base_raw}")

classify_raw=$(sbatch \
  --parsable \
  --chdir="${repo_root}" \
  --job-name=inflect-modality-classify \
  --dependency="afterok:${pilot_id}:${base_id}" \
  --time=02:00:00 \
  --mem=16G \
  --export=ALL,INFLECT_MODALITY_MODE=classify \
  "${job_script}")
classify_id=$(job_id "${classify_raw}")

refine_raw=$(sbatch \
  --parsable \
  --chdir="${repo_root}" \
  --job-name=inflect-modality-refine \
  --dependency="afterok:${classify_id}" \
  --array="0-$((array_count - 1))%${array_limit}" \
  --time=24:00:00 \
  --mem=8G \
  --export="ALL,INFLECT_MODALITY_MODE=refine,INFLECT_ARRAY_COUNT=${array_count}" \
  "${job_script}")
refine_id=$(job_id "${refine_raw}")

summarize_raw=$(sbatch \
  --parsable \
  --chdir="${repo_root}" \
  --job-name=inflect-modality-summarize \
  --dependency="afterok:${refine_id}" \
  --time=06:00:00 \
  --mem=32G \
  --export=ALL,INFLECT_MODALITY_MODE=summarize \
  "${job_script}")
summarize_id=$(job_id "${summarize_raw}")

manifest="${work_dir}/submitted-jobs.tsv"
{
  printf 'stage\t%s\n' "${stage_id}"
  printf 'pilot\t%s\n' "${pilot_id}"
  printf 'base\t%s\n' "${base_id}"
  printf 'classify\t%s\n' "${classify_id}"
  printf 'refine\t%s\n' "${refine_id}"
  printf 'summarize\t%s\n' "${summarize_id}"
} > "${manifest}"

printf 'Submitted checkpointed modality workflow:\n'
sed 's/^/  /' "${manifest}"
