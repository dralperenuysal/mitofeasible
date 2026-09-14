#!/bin/bash
#SBATCH --job-name=mtcov_indel
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=03:00:00
#SBATCH --output=/arf/scratch/suysal/mtcovmap/logs/indel_%j.out
#SBATCH --error=/arf/scratch/suysal/mtcovmap/logs/indel_%j.err
set -euo pipefail
ROOT=${MTCOV_ROOT:-/arf/scratch/suysal/mtcovmap}
REPO=${MTCOV_REPO:-/arf/home/suysal/mtcovmap}
mkdir -p "$ROOT/logs"
apptainer exec --bind /arf "$ROOT/mtcovmap.sif" python3 \
    "$REPO/scripts/11_indel_at_homopolymers.py" \
    --config "$REPO/config/params.yaml" \
    --samples "$REPO/config/samples.tsv" \
    --out "$REPO/results/tables"
