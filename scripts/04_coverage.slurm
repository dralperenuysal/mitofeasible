#!/bin/bash
# Phase 4 - per-base coverage, one array task per sample.
#
#   cd /arf/scratch/suysal/mtcovmap
#   sbatch --array=1-60%10 ~/mtcovmap/scripts/04_coverage.slurm
#
#SBATCH --job-name=mtcov_cov
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=04:00:00
#SBATCH --output=/arf/scratch/suysal/mtcovmap/logs/cov_%A_%a.out
#SBATCH --error=/arf/scratch/suysal/mtcovmap/logs/cov_%A_%a.err

set -euo pipefail
ROOT=${MTCOV_ROOT:-/arf/scratch/suysal/mtcovmap}   # scratch: big intermediates
REPO=${MTCOV_REPO:-/arf/home/suysal/mtcovmap}       # this checkout, on the cluster
APPT="apptainer exec --bind /arf $ROOT/mtcovmap.sif"
ROW=${SLURM_ARRAY_TASK_ID:?run under sbatch --array}
mkdir -p "$ROOT/logs" "$ROOT/coverage"
cd "$REPO"

RUN=$($APPT python -c "
import pandas as pd; print(pd.read_csv('config/samples.tsv', sep='\t',
dtype=str).iloc[$ROW - 1]['run_accession'])")
[ -f "$ROOT/coverage/${RUN}_per_base.parquet" ] && { echo "have $RUN"; exit 0; }

# Same lock idiom as Phase 3, so the array can be twinned across partitions.
LOCK=$ROOT/coverage/.lock_$RUN
mkdir "$LOCK" 2>/dev/null || { echo "$RUN taken by a twin job"; exit 0; }
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

$APPT python scripts/04_coverage.py --config config/params.yaml \
    --samples config/samples.tsv --row "$ROW"
