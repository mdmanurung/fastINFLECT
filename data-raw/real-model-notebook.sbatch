#!/usr/bin/env bash
#SBATCH --cpus-per-task=2
#SBATCH --output=data-raw/.real-model-notebook-work/notebook-%j.out
#SBATCH --error=data-raw/.real-model-notebook-work/notebook-%j.err

set -euo pipefail

repo_root=${SLURM_SUBMIT_DIR:?Submit this job from the fastINFLECT repository}
jupyter_bin=${INFLECT_NOTEBOOK_JUPYTER:-/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/codex/bin/jupyter}
rscript_bin=${INFLECT_NOTEBOOK_RSCRIPT:-/exports/archive/hg-funcgenom-research/mdmanurung/conda/envs/R4_51/bin/Rscript}
work_dir=${INFLECT_NOTEBOOK_WORKDIR:-${repo_root}/data-raw/.real-model-notebook-work}
notebook=${repo_root}/data-raw/test_inflect.ipynb
executed_name=test_inflect.executed.ipynb
executed=${work_dir}/${executed_name}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_LIBS=${repo_root}/data-raw/.real-model-lib:/exports/para-lipg-hpc/mdmanurung/R/4.5

mkdir -p "${work_dir}"
cd "${repo_root}"

"${rscript_bin}" data-raw/summarize-bmv-k-selection.R

"${jupyter_bin}" nbconvert \
  --to notebook \
  --execute \
  --ExecutePreprocessor.kernel_name=ir \
  --ExecutePreprocessor.timeout=-1 \
  --output-dir="${work_dir}" \
  --output="${executed_name}" \
  "${notebook}"

jq -e '
  [.cells[] | select(.cell_type == "code") | .execution_count] |
  length > 0 and all(. != null)
' "${executed}" >/dev/null

mv "${executed}" "${notebook}"
printf 'Executed notebook: %s\n' "${notebook}"
