#!/bin/bash
# Phase 2 retrieval, run on the ARF login node: the work is network I/O, not CPU,
# so a 40-core compute allocation would idle.
# Restartable: verified files are skipped, partial ones resume (curl -C -).
# Stream count is a tuning knob, not a constant: ENA throttles per connection,
# so the useful number is whatever measurement says, not a default.
set -uo pipefail
cd "$(dirname "$0")/.."
source /etc/profile.d/modules.sh 2>/dev/null || true
module load miniconda3
n=$(($(wc -l < config/samples.tsv) - 1))
seq 1 "$n" | xargs -P "${1:-4}" -I{} sh -c \
  'python scripts/02_download.py --download --row {} >> logs/dl_{}.log 2>&1 || echo "FAILED row {}"'
echo "all rows attempted"
