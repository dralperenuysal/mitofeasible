#!/bin/bash
# Chemistry-control cohort (PRJNA1086804, bulk DLPFC, poly(A) vs rRNA-depleted).
# Downloads run on the login node: network I/O is a file operation, not compute
# (truba.md). Resumable and checksum-verified, so an interrupted run is re-runnable.
#
#   nohup bash ~/mtcovmap/scripts/10_fetch_chem.sh > ~/mtcovmap/logs/fetch_chem.log 2>&1 &
set -uo pipefail
: "${MTCOV_REPO:?set MTCOV_REPO to this checkout (run_all.sh exports it)}"
source "$MTCOV_REPO/scripts/_env.sh"
DEST=$ROOT/fastq_chem
mkdir -p "$DEST" "$REPO/logs"

# Single-instance lock. Launching this over a flaky ssh can start a second copy
# while the first is still alive; two curls writing one path corrupt it, and the
# md5 check then fails the file forever. Atomic mkdir, same idiom as the aligner.
LOCK=$DEST/.fetch.lock
if ! mkdir "$LOCK" 2>/dev/null; then
    echo "another fetch holds $LOCK -- exiting (remove it if no fetch is running)"
    exit 0
fi
# Only the main shell may release the lock. Background `fetch` jobs run in
# subshells that inherit this trap, so an unguarded one is removed by the first
# download that finishes - which is how three copies ended up on one file.
trap '[ "$BASHPID" = "$$" ] && rmdir "$LOCK" 2>/dev/null' EXIT

# Three streams, not four: ENA refuses connections on port 21 under a heavier burst.
STREAMS=3

fetch() {                       # url  md5  destination
    local url=$1 want=$2 f=$3
    if [ -s "$f" ] && [ "$(md5sum "$f" | cut -d' ' -f1)" = "$want" ]; then
        echo "ok (cached)  $(basename "$f")"; return 0
    fi
    for resume in 1 0; do
        local args=(-fsSL --retry 10 --retry-delay 15 --retry-all-errors
                    --speed-time 60 --speed-limit 10000 -o "$f")
        [ "$resume" = 1 ] || rm -f "$f"
        [ "$resume" = 1 ] && args=(-C - "${args[@]}")
        curl "${args[@]}" "ftp://$url" && \
          [ "$(md5sum "$f" | cut -d' ' -f1)" = "$want" ] && {
            echo "ok           $(basename "$f")"; return 0; }
        echo "retry (full) $(basename "$f")"
    done
    echo "FAILED       $(basename "$f")"; return 1
}
export -f fetch

# One line per file, so xargs can run a few streams in parallel.
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next}
            {n=split($h["fastq_ftp"],u,";"); split($h["fastq_md5"],m,";");
             for(i=1;i<=n;i++) print u[i]"\t"m[i]}' "$REPO/config/samples_chem.tsv" \
| while IFS=$'\t' read -r url md5; do
      printf '%s\t%s\t%s\n' "$url" "$md5" "$DEST/$(basename "$url")"
  done > "$DEST/.manifest"

wc -l < "$DEST/.manifest" | xargs echo "files to fetch:"
# A plain read loop, four streams at a time. (xargs -I with a
# tab-bearing replacement string mangled the fields.)
while IFS=$'\t' read -r url md5 dest; do
    while [ "$(jobs -rp | wc -l)" -ge "$STREAMS" ]; do wait -n; done
    fetch "$url" "$md5" "$dest" &
done < "$DEST/.manifest"
wait
echo "=== download stage done: $(date) ==="
ls -1 "$DEST"/*.fastq.gz 2>/dev/null | wc -l | xargs echo "fastq files present:"
