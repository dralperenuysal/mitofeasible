#!/bin/bash
#SBATCH --job-name=mtcov_indel
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=03:00:00
#SBATCH --output=logs/indel_%j.out
#SBATCH --error=logs/indel_%j.err
set -euo pipefail
: "${MTCOV_REPO:?set MTCOV_REPO to this checkout (run_all.sh exports it)}"
source "$MTCOV_REPO/scripts/_env.sh"
mkdir -p "$ROOT/logs"
pyrun "$REPO/scripts/11_indel_at_homopolymers.py" \
    --config "$REPO/config/params.yaml" \
    --samples "$REPO/config/samples.tsv" \
    --out "$REPO/results/tables"
