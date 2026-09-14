#!/bin/bash
# Phase 3 - align one sample twice (A: NUMTs intact, B: NUMTs masked), then
# realign its chrM reads to the shifted chrM so the control region is not cut by
# the linear reference's seam (scripts/03_circularity_decision.md).
#
#   cd /arf/scratch/suysal/mtcovmap
#   sbatch --array=1-60%6 ~/mtcovmap/scripts/03_align.sh
#
#SBATCH --job-name=mtcov_aln
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=12:00:00
#SBATCH --output=/arf/scratch/suysal/mtcovmap/logs/aln_%A_%a.out
#SBATCH --error=/arf/scratch/suysal/mtcovmap/logs/aln_%A_%a.err

set -euo pipefail
# STAR's coordinate sort opens outBAMsortingThreadN x outBAMsortingBinsN temp
# files at once. With the default (one sorting thread per aligner thread, 50 bins)
# that is 1000 files at 20 threads and 2800 at 56 - at or over the open-file
# limit, and it fails only after the genome is loaded. Raise the soft limit and
# cap the sorting threads; both are done here so A, B and the shifted pass share
# the setting.
ulimit -n "$(ulimit -Hn)" 2>/dev/null || true
ROOT=${MTCOV_ROOT:-/arf/scratch/suysal/mtcovmap}   # scratch: big intermediates
REPO=${MTCOV_REPO:-/arf/home/suysal/mtcovmap}       # this checkout, on the cluster
# A second cohort (the chemistry control) has its own manifest, its own flat
# FASTQ directory and needs only the unmasked arm. Defaults reproduce the
# original behaviour exactly, so the main cohort is unaffected.
MANIFEST=${MTCOV_MANIFEST:-$REPO/config/samples.tsv}
ARMS=${MTCOV_ARMS:-A B}
SIF=$ROOT/mtcovmap.sif
APPT="apptainer exec --bind /arf $SIF"
ROW=${SLURM_ARRAY_TASK_ID:?set SLURM_ARRAY_TASK_ID or run under sbatch --array}
mkdir -p "$ROOT/logs" "$ROOT/bam" "$ROOT/stats"

# One line of the manifest. Reading it here keeps sample identity out of the
# script and in config/samples.tsv, where Phase 2 put it.
# The heredoc is quoted, so the checkout path is passed as an argument rather
# than interpolated - the script must not depend on the shell expanding inside it.
read -r RUN LAYOUT THREADS MULTIMAP MISMATCH SJMIN SORTRAM EXTRA < <($APPT python - "$ROW" "$REPO" "$MANIFEST" <<'PY'
import sys, pandas as pd, yaml
repo = sys.argv[2]
row = pd.read_csv(sys.argv[3], sep="\t", dtype=str).iloc[int(sys.argv[1]) - 1]
c = yaml.safe_load(open(f"{repo}/config/params.yaml"))["alignment"]
print(row["run_accession"], row["library_layout"], c["threads"],
      c["out_filter_multimap_nmax"], c["out_filter_mismatch_nover_lmax"],
      c["align_sjdb_overhang_min"], c["limit_bam_sort_ram"], c["star_extra"] or "-")
PY
)
[ "$EXTRA" = "-" ] && EXTRA=""
# config/params.yaml carries barbun's minimum; a twin job on a wider partition
# gets a bigger allocation, and leaving it idle would be the real waste. STAR's
# alignments do not depend on the thread count, only the wall time does.
THREADS=${SLURM_CPUS_ON_NODE:-$THREADS}
# Atomic lock per sample, so the same row can be submitted to two partitions at
# once (truba.md: whichever queue opens first wins) without two jobs writing the
# same BAM. The loser exits without touching anything.
LOCK=$ROOT/bam/.lock_$RUN
if ! mkdir "$LOCK" 2>/dev/null; then
    echo "$RUN is already being aligned by a twin job -- exiting."
    exit 0
fi
trap 'rmdir "$LOCK" 2>/dev/null' EXIT

# The main cohort keeps one directory per run; the chemistry cohort is flat.
FQ=${MTCOV_FQDIR:-$ROOT/fastq/$RUN}
if [ "$LAYOUT" = "PAIRED" ]; then
    READS="$FQ/${RUN}_1.fastq.gz $FQ/${RUN}_2.fastq.gz"
else
    READS="$FQ/${RUN}.fastq.gz"
fi

# Identical for A and B by construction: one variable, used twice.
SORT_THREADS=$(( THREADS < 8 ? THREADS : 8 ))
star_args=(--runThreadN "$THREADS" --outBAMsortingThreadN "$SORT_THREADS"
           --outFilterMultimapNmax "$MULTIMAP"
           --outFilterMismatchNoverLmax "$MISMATCH"
           --alignSJDBoverhangMin "$SJMIN"
           --limitBAMsortRAM "$SORTRAM"
           --outSAMtype BAM SortedByCoordinate
           --outSAMattributes NH HI AS nM NM MD)
# readFilesCommand is added per call, not here: both passes read .gz, but the
# shifted pass reads files this script writes, so the coupling stays visible.

for g in $ARMS; do
    mt=$ROOT/bam/${RUN}_${g}_chrM.bam
    spre=$ROOT/bam/${RUN}_${g}_shifted_
    # The genome-wide BAM is deleted at the end of each arm, so resume on the
    # chrM products - the only outputs that survive.
    [ -f "$mt.bai" ] && [ -f "${spre}Aligned.sortedByCoord.out.bam.bai" ] && {
        echo "have $RUN $g"; continue; }

    pre=$ROOT/bam/${RUN}_${g}_
    bam=${pre}Aligned.sortedByCoord.out.bam
    $APPT STAR --genomeDir "$ROOT/index/genome_$g" \
        --readFilesIn $READS --outFileNamePrefix "$pre" \
        "${star_args[@]}" --readFilesCommand zcat $EXTRA
    $APPT samtools index -@ "$THREADS" "$bam"

    # Duplicate rate is an alignment-level number AGENTS.md asks for. Marked,
    # never removed: PCR duplicates and genuinely deep chrM coverage are not
    # distinguishable here, and dropping them would bias depth downwards exactly
    # where coverage is highest.
    $APPT samtools view -b -@ "$THREADS" "$bam" chrM \
        | $APPT samtools collate -u -O -@ "$THREADS" - \
        | $APPT samtools fixmate -m -u -@ "$THREADS" - - \
        | $APPT samtools sort -u -@ "$THREADS" - \
        | $APPT samtools markdup -@ "$THREADS" \
              -f "$ROOT/stats/${RUN}_${g}.markdup" - "$mt"
    $APPT samtools index "$mt"

    # Alignment-level numbers, recorded before anything is filtered.
    {
        printf '%s\t%s\ttotal\t%s\n' "$RUN" "$g" \
            "$($APPT samtools view -c -@ "$THREADS" "$bam")"
        printf '%s\t%s\tchrM\t%s\n' "$RUN" "$g" \
            "$($APPT samtools view -c -@ "$THREADS" "$mt")"
    } > "$ROOT/stats/${RUN}_${g}.counts"
    $APPT samtools flagstat -@ "$THREADS" "$mt" > "$ROOT/stats/${RUN}_${g}.flagstat"

    # Shifted pass: realign the chrM reads only. The seam moves into MT-CO2, so
    # the control region is read off a continuous stretch.
    fq=$ROOT/bam/${RUN}_${g}_chrM
    if [ "$LAYOUT" = "PAIRED" ]; then
        $APPT samtools collate -u -O -@ "$THREADS" "$mt" \
            | $APPT samtools fastq -n -@ "$THREADS" \
                  -1 "${fq}_1.fq.gz" -2 "${fq}_2.fq.gz" -0 /dev/null -s /dev/null -
        SHIFT_READS="${fq}_1.fq.gz ${fq}_2.fq.gz"
    else
        $APPT samtools fastq -n -@ "$THREADS" "$mt" | gzip > "${fq}.fq.gz"
        SHIFT_READS="${fq}.fq.gz"
    fi
    $APPT STAR --genomeDir "$ROOT/index/chrM_shifted" \
        --readFilesIn $SHIFT_READS --outFileNamePrefix "$spre" \
        "${star_args[@]}" --readFilesCommand zcat $EXTRA
    $APPT samtools index "${spre}Aligned.sortedByCoord.out.bam"
    rm -f "${fq}"_*.fq.gz "${fq}.fq.gz" "$bam" "$bam.bai"
done

echo "done $RUN"
