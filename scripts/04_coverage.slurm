#!/bin/bash
# Phase 4 - per-base coverage, one array task per sample.
#
#   cd "$MTCOV_ROOT"
#   sbatch --array=1-60%10 ~/mtcovmap/scripts/04_coverage.slurm
#
#SBATCH --job-name=mtcov_cov
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=04:00:00
#SBATCH --output=logs/cov_%A_%a.out
#SBATCH --error=logs/cov_%A_%a.err

set -euo pipefail
: "${MTCOV_REPO:?set MTCOV_REPO to this checkout (run_all.sh exports it)}"
source "$MTCOV_REPO/scripts/_env.sh"
ROW=${SLURM_ARRAY_TASK_ID:?run under sbatch --array}
mkdir -p "$ROOT/logs" "$ROOT/coverage"
cd "$REPO"

RUN=$($APPT python3 -c "
import pandas as pd; print(pd.read_csv('config/samples.tsv', sep='\t',
dtype=str).iloc[$ROW - 1]['run_accession'])")
[ -f "$ROOT/coverage/${RUN}_per_base.parquet" ] && { echo "have $RUN"; exit 0; }

# Same lock idiom as Phase 3, so the array can be twinned across partitions.
LOCK=$ROOT/coverage/.lock_$RUN
mkdir "$LOCK" 2>/dev/null || { echo "$RUN taken by a twin job"; exit 0; }
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

$APPT python3 scripts/04_coverage.py --config config/params.yaml \
    --samples config/samples.tsv --row "$ROW"
