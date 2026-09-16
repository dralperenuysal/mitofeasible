#!/bin/bash
# Coverage + allele counts for one chemistry-cohort sample. Run as an array so
# the 38 samples go in parallel: the same work in a sequential shell loop took
# ~70 min for 24 samples, almost all of it in the poly(A) libraries, whose chrM
# read counts are ~18x the rRNA-depleted ones.
#
#   cd "$MTCOV_ROOT"
#   sbatch --array=1-38 scripts/15_chem_counts.slurm
#
#SBATCH --job-name=mtcov_chemcnt
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=02:00:00
#SBATCH --output=logs/chemcnt_%A_%a.out
#SBATCH --error=logs/chemcnt_%A_%a.err
set -euo pipefail
: "${MTCOV_REPO:?set MTCOV_REPO to this checkout (run_all.sh exports it)}"
source "$MTCOV_REPO/scripts/_env.sh"
MANIFEST=$REPO/config/samples_chem.tsv
ROW=${SLURM_ARRAY_TASK_ID:?}

RUN=$($APPT python3 -c "
import pandas as pd, sys
print(pd.read_csv('$MANIFEST', sep='\t', dtype=str).iloc[$ROW - 1]['run_accession'])")

if [ ! -f "$ROOT/coverage/${RUN}_per_gene.tsv" ]; then
    $APPT python3 "$REPO/scripts/04_coverage.py" --config "$REPO/config/params.yaml" \
        --samples "$MANIFEST" --row "$ROW" --stranded false --arms A
else
    echo "have coverage $RUN"
fi

if [ ! -f "$ROOT/alleles/${RUN}_alleles.parquet" ]; then
    $APPT python3 "$REPO/scripts/07a_allele_counts.py" --config "$REPO/config/params.yaml" \
        --samples "$MANIFEST" --row "$ROW" --arms A
else
    echo "have alleles $RUN"
fi
echo "done $RUN"
