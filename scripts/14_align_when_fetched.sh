#!/bin/bash
# Wait for the chemistry cohort download to finish, then align whatever is
# complete and not yet aligned. Runs unattended on the login node so the queue
# time is not lost while nobody is watching.
#
#   setsid nohup bash ~/mtcovmap/scripts/14_align_when_fetched.sh \
#       > ~/mtcovmap/logs/align_when_fetched.log 2>&1 < /dev/null &
set -uo pipefail
: "${MTCOV_REPO:?set MTCOV_REPO to this checkout (run_all.sh exports it)}"
source "$MTCOV_REPO/scripts/_env.sh"
LOG=${1:-$REPO/logs/fetch_chem2.log}
MANIFEST=$REPO/config/samples_chem.tsv
FQ=$ROOT/fastq_chem
PARTITION=${MTCOV_PARTITION:-barbun}     # production queue, not debug

# Single instance, same idiom (and same subshell caveat) as the fetch script.
LOCK=$ROOT/.align_watch.lock
mkdir "$LOCK" 2>/dev/null || { echo "another watcher holds $LOCK -- exiting"; exit 0; }
trap '[ "$BASHPID" = "$$" ] && rmdir "$LOCK" 2>/dev/null' EXIT

echo "watching $LOG  ($(date))"
for _ in $(seq 1 480); do                # 480 x 60s = 8h ceiling
    grep -q "download stage done" "$LOG" 2>/dev/null && break
    sleep 60
done
if ! grep -q "download stage done" "$LOG" 2>/dev/null; then
    echo "download never reported done; not submitting"; exit 1
fi

failed=$(grep -c "^FAILED" "$LOG" 2>/dev/null || echo 0)
echo "download finished with $failed failed file(s)"

# Rows whose reads are both present and whose chrM BAM does not exist yet. Done
# in awk against the manifest so the row numbers match what 03_align.sh expects.
rows=$(awk -F'\t' -v fq="$FQ" -v bam="$ROOT/bam" '
    NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next }
    {
        r=$h["run_accession"]; n=NR-1
        cmd="test -s " fq "/" r "_1.fastq.gz && test -s " fq "/" r "_2.fastq.gz"
        if (system(cmd)!=0) next
        if (system("test -f " bam "/" r "_A_chrM.bam.bai")==0) next
        printf "%s%d", (out++ ? "," : ""), n
    }' "$MANIFEST")

if [ -z "$rows" ]; then
    echo "nothing left to align"; exit 0
fi
echo "submitting rows: $rows"
cd "$ROOT" || exit 1                     # sbatch is refused outside /arf/scratch
sbatch --array="$rows" -p "$PARTITION" --time=04:00:00 \
    --export=ALL,MTCOV_MANIFEST=$MANIFEST,MTCOV_FQDIR=$FQ,MTCOV_ARMS=A \
    "$REPO/scripts/03_align.sh"
echo "done ($(date))"
