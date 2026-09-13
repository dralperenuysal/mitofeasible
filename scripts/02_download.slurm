#!/bin/bash
# Phase 2 retrieval on TRUBA. One array task per sample; each verifies the ENA
# md5 and is safe to re-run (curl -C - resumes, verified files are skipped).
#
#   sbatch --array=1-51%4 scripts/02_download.slurm
#
# %4 caps concurrency: ENA throttles, and four streams already saturate.
#SBATCH --job-name=mtcov_dl
#SBATCH --output=logs/dl_%A_%a.out
#SBATCH --error=logs/dl_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=08:00:00
##SBATCH --partition=<set this>   # fill in from `sinfo`; TRUBA queue names vary

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mitofeasible
python scripts/02_download.py --download --row "${SLURM_ARRAY_TASK_ID}"
