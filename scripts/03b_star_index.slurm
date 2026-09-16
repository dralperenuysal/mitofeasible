#!/bin/bash
# Phase 3b - build the STAR indices. One array task per index, so A, B and the
# shifted chrM build concurrently instead of end to end; A and B are the long
# poles (~1 h each) and nothing couples them.
#
#   cd "$MTCOV_ROOT"
#   sbatch --array=0-2 ~/mtcovmap/scripts/03b_star_index.slurm
#
# Task 0 also prepares the references, and 1 and 2 wait for it: all three read
# ref/, and only one may write it.
#SBATCH --job-name=mtcov_idx
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=08:00:00
#SBATCH --output=logs/idx_%A_%a.out
#SBATCH --error=logs/idx_%A_%a.err

set -euo pipefail
: "${MTCOV_REPO:?set MTCOV_REPO to this checkout (run_all.sh exports it)}"
source "$MTCOV_REPO/scripts/_env.sh"
TASK=${SLURM_ARRAY_TASK_ID:-0}
mkdir -p "$ROOT/logs" "$ROOT/index"

TARGETS=(genome_A genome_B chrM_shifted)
NAME=${TARGETS[$TASK]}
OUT=$ROOT/index/$NAME
READY=$ROOT/ref/.references_ready

if [ "$TASK" -eq 0 ]; then
    $APPT python3 "$REPO/scripts/03a_prepare_reference.py" \
        --config "$REPO/config/params.yaml" --ref-dir "$ROOT/ref"
    touch "$READY"
else
    # References are written once, by task 0. Wait rather than race.
    for _ in $(seq 1 240); do [ -f "$READY" ] && break; sleep 30; done
    [ -f "$READY" ] || { echo "references never appeared"; exit 1; }
fi

# An atomic mkdir lock, so a resubmitted array cannot have two tasks building the
# same index into the same directory (truba.md).
LOCK=$ROOT/index/.lock_$NAME
if ! mkdir "$LOCK" 2>/dev/null; then
    echo "$NAME is being built by another task -- exiting."
    exit 0
fi
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

# -s, not -f: an empty SA is a placeholder or a truncated build, not an index.
if [ -s "$OUT/SA" ]; then
    echo "have $NAME"
    exit 0
fi
# Build into a staging directory and move it into place at the end, so a
# half-built index is never visible under its final name - a killed job would
# otherwise leave a directory that looks finished to the next reader.
BUILD=$ROOT/index/.building_$NAME
rm -rf "$BUILD"; mkdir -p "$BUILD"

read -r THREADS OVERHANG < <($APPT python3 -c "
import yaml; c = yaml.safe_load(open('$REPO/config/params.yaml'))['alignment']
print(c['threads'], c['sjdb_overhang'])")

if [ "$NAME" = "chrM_shifted" ]; then
    # 16.6 kb: STAR's suffix-array index must be shrunk or it builds a broken one.
    $APPT STAR --runMode genomeGenerate --runThreadN "$THREADS" \
        --genomeDir "$BUILD" --genomeFastaFiles "$ROOT/ref/chrM_shifted.fa" \
        --genomeSAindexNbases 6 --outFileNamePrefix "$BUILD/"
else
    $APPT STAR --runMode genomeGenerate --runThreadN "$THREADS" \
        --genomeDir "$BUILD" --genomeFastaFiles "$ROOT/ref/${NAME}.fa" \
        --sjdbGTFfile "$ROOT/ref/annotation.gtf" \
        --sjdbOverhang "$OVERHANG" --outFileNamePrefix "$BUILD/"
fi

rm -rf "$OUT"
mv "$BUILD" "$OUT"
echo "index ready: $NAME"
