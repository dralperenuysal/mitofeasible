#!/bin/bash
# Phase 6 - haplogroup profiles, one array task per sample.
#
#   cd /arf/scratch/suysal/mtcovmap
#   sbatch --array=1-60%10 ~/mtcovmap/scripts/06_haplogroup.slurm
#
#SBATCH --job-name=mtcov_alc
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --time=03:50:00
#SBATCH --output=/arf/scratch/suysal/mtcovmap/logs/alc_%A_%a.out
#SBATCH --error=/arf/scratch/suysal/mtcovmap/logs/alc_%A_%a.err

set -euo pipefail
ROOT=/arf/scratch/suysal/mtcovmap
REPO=/arf/home/suysal/mtcovmap
APPT="apptainer exec --bind /arf $ROOT/mtcovmap.sif"
ROW=${SLURM_ARRAY_TASK_ID:?run under sbatch --array}
mkdir -p "$ROOT/logs" "$ROOT/alleles"
cd "$REPO"

RUN=$($APPT python -c "
import pandas as pd; print(pd.read_csv('config/samples.tsv', sep='\t',
dtype=str).iloc[$ROW - 1]['run_accession'])")
[ -f "$ROOT/alleles/${RUN}_alleles.parquet" ] && { echo "have $RUN"; exit 0; }

LOCK=$ROOT/alleles/.lock_$RUN
mkdir "$LOCK" 2>/dev/null || { echo "$RUN taken by a twin job"; exit 0; }
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

$APPT python scripts/07a_allele_counts.py --row "$ROW"
