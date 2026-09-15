#!/bin/bash
# The order the phases actually ran in, with the exact invocations.
#
# The README says what each phase does; this says how it was called, which is
# the part that was previously only in a terminal history. Several invocations
# are not guessable from the scripts alone - the chemistry cohort in particular
# needs three environment overrides and two flags.
#
#   ./run_all.sh              # print the plan and exit (default; runs nothing)
#   ./run_all.sh 3 4 5        # run those phases
#   ./run_all.sh all          # everything that does not need human review
#
# Phases 0 and 1 are deliberately absent from `all`: cohort selection is a
# scientific decision, and both halt for review (README, Pipeline).
set -uo pipefail

ROOT=${MTCOV_ROOT:-/arf/scratch/suysal/mtcovmap}   # scratch: big intermediates
REPO=${MTCOV_REPO:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
SIF=${MTCOV_SIF:-$ROOT/mtcovmap.sif}
PART=${MTCOV_PARTITION:-barbun}                    # production queue; NOT debug
# Everything heavy runs inside the container. The selftests do not: they need no
# data and no bioinformatics tools, so on a laptop with the Python dependencies
# installed they should still run. Fall back to the interpreter when there is no
# apptainer rather than failing in a way that looks like a broken checkout.
if command -v apptainer >/dev/null 2>&1; then
    APPT="apptainer exec --bind /arf $SIF"
else
    APPT=""
    echo "note: apptainer not found; running scripts directly (selftests only)" >&2
fi
N_MAIN=60                                          # 30 LCL + 30 muscle
N_CHEM=38                                          # 19 blocks x 2 chemistries

say() { printf '\n\033[1m== %s ==\033[0m\n' "$*"; }
run() { echo "+ $*"; [ "${DRY:-0}" = 1 ] || "$@"; }
py()  { if [ -n "$APPT" ]; then $APPT python3 "$@"; else python3 "$@"; fi; }

# sbatch is refused outside /arf/scratch on this cluster (truba.md), so phases
# that submit jobs chdir there first. Done lazily: `plan` and `selftest` must
# work from a plain clone with no cluster attached.
need_root() {
    cd "$ROOT" 2>/dev/null || {
        echo "MTCOV_ROOT=$ROOT not reachable -- set it, or use DRY=1" >&2
        [ "${DRY:-0}" = 1 ] || exit 1; }
}

phase_env() {
    say "env: build the analysis container"
    echo "  On TRUBA:   bash $REPO/scripts/truba_build_container.sh"
    echo "  Elsewhere:  apptainer build mitofeasible.sif $REPO/env/mitofeasible.def"
    echo "  Versions that produced the published numbers: $REPO/env/versions.lock"
}

phase_2() {
    need_root
    say "2  manifest, then download"
    # Without --download this writes config/samples.tsv and stops, which is the
    # point: the manifest is worth reading before 200 GB moves.
    run $APPT python3 "$REPO/scripts/02_download.py" --config "$REPO/config/params.yaml"
    echo "  review config/samples.tsv, then fetch (%4 caps concurrency; ENA throttles):"
    run sbatch --array=1-$N_MAIN%4 "$REPO/scripts/02_download.slurm"
}

phase_3() {
    need_root
    say "3  reference, index, then align every sample twice (A intact, B masked)"
    run $APPT python3 "$REPO/scripts/03a_prepare_reference.py" --config "$REPO/config/params.yaml"
    run sbatch --array=0-2 "$REPO/scripts/03b_star_index.slurm"
    echo "  wait for the three indexes, then:"
    run sbatch --array=1-$N_MAIN%6 -p "$PART" "$REPO/scripts/03_align.sh"
}

phase_4() { need_root; say "4  coverage";       run sbatch --array=1-$N_MAIN%10 -p "$PART" "$REPO/scripts/04_coverage.slurm"; }
phase_5() { need_root; say "5  NUMT delta";     run $APPT python3 "$REPO/scripts/05_numt_delta.py" --config "$REPO/config/params.yaml"; }
phase_6() {
    need_root
    say "6  haplogroups, then compare against the donors' own DNA"
    run sbatch --array=1-$N_MAIN%10 -p "$PART" "$REPO/scripts/06_haplogroup.slurm"
    run $APPT python3 "$REPO/scripts/06_haplogroup.py" classify --config "$REPO/config/params.yaml"
    run $APPT python3 "$REPO/scripts/06b_truth_compare.py" --config "$REPO/config/params.yaml"
}
phase_7() {
    need_root
    say "7  allele counts, then the feasibility map"
    run sbatch --array=1-$N_MAIN%10 -p "$PART" "$REPO/scripts/07a_allele_counts.slurm"
    run $APPT python3 "$REPO/scripts/07_feasibility_map.py" --config "$REPO/config/params.yaml"
}
phase_8() { need_root; say "8  does the map transfer?"; run $APPT python3 "$REPO/scripts/08_replication.py" --config "$REPO/config/params.yaml"; }

phase_11() {
    need_root
    say "11  is the poly-C noise length ambiguity? (reads existing BAMs)"
    run sbatch -p "$PART" "$REPO/scripts/11_indel.slurm"
}

phase_chem() {
    need_root
    say "12  chemistry cohort: poly(A) vs rRNA depletion"
    # This cohort is not in config/params.yaml on purpose: phase 8 picks the
    # replication cohort as "the study that is not the primary one", and a third
    # entry there would silently change that. It is carried by its own manifest
    # and three overrides instead.
    local M=$REPO/config/samples_chem.tsv
    local FQ=$ROOT/fastq_chem
    echo "  fetch (login node, resumable, md5-verified):"
    run bash "$REPO/scripts/10_fetch_chem.sh"
    echo "  align: unmasked arm only - the NUMT contrast is not what this cohort answers"
    run sbatch --array=1-$N_CHEM -p "$PART" \
        --export=ALL,MTCOV_MANIFEST=$M,MTCOV_FQDIR=$FQ,MTCOV_ARMS=A \
        "$REPO/scripts/03_align.sh"
    echo "  coverage + allele counts, one array (sequential loop took ~70 min for 24):"
    run sbatch --array=1-$N_CHEM -p "$PART" "$REPO/scripts/15_chem_counts.slurm"
    echo "  paired analysis over the 19 tissue blocks:"
    run $APPT python3 "$REPO/scripts/12_chemistry.py" \
        --config "$REPO/config/params.yaml" --samples "$M"
}

phase_9() {
    need_root
    say "9  figures"
    run $APPT python3 "$REPO/scripts/09_figures.py" --config "$REPO/config/params.yaml"
    run $APPT python3 "$REPO/scripts/09b_circularity.py" --config "$REPO/config/params.yaml"
}

selftest() {
    say "selftests (no data needed; run these first on a fresh clone)"
    local fail=0
    for f in "$REPO"/scripts/*.py; do
        grep -q selftest "$f" || continue
        case $(basename "$f") in
            # $APPT is a command plus flags, so it must word-split: quoting it
            # makes the whole string one command name.
            06_haplogroup.py) py "$f" selftest  >/dev/null 2>&1 || { echo "  FAIL $(basename "$f")"; fail=1; } ;;
            *)                py "$f" --selftest >/dev/null 2>&1 || { echo "  FAIL $(basename "$f")"; fail=1; } ;;
        esac
    done
    [ "$fail" = 0 ] && echo "  all selftests pass" || echo "  SOME SELFTESTS FAILED" >&2
    return $fail
}

plan() {
    cat <<'TXT'
Phases, in order. Phases 0 and 1 halt for human review and are not automated.

  env    build the container            (TRUBA script, or env/mitofeasible.def)
  2      download FASTQ                 ~200 GB, md5-verified against ENA
  3      reference, index, align        60 samples x 2 arms + a shifted chrM pass
  4      coverage                       per-base and per-gene
  5      NUMT delta                     where masking changes the picture
  6      haplogroups                    checked against the donors' own DNA
  7      allele counts, feasibility map the central result
  8      replication                    does the map transfer?
  11     indel check at the poly-C tracts
  12     chemistry cohort               separate manifest + 3 overrides
  9      figures

  ./run_all.sh selftest      verify the logic without any data
  ./run_all.sh 2 3 4         run specific phases
  ./run_all.sh all           2,3,4,5,6,7,8,11,chem,9
  DRY=1 ./run_all.sh all     print every command without running it

Paths come from MTCOV_ROOT, MTCOV_REPO, MTCOV_SIF, MTCOV_PARTITION.
Array jobs are submitted, not waited on: check `squeue` before the next phase.
TXT
}

[ $# -eq 0 ] && { plan; exit 0; }
for p in "$@"; do
    case $p in
        all)      for q in 2 3 4 5 6 7 8 11 chem 9; do "phase_$q"; done ;;
        selftest) selftest ;;
        plan)     plan ;;
        env|2|3|4|5|6|7|8|9|11|chem) "phase_$p" ;;
        *) echo "unknown phase: $p" >&2; plan; exit 2 ;;
    esac
done
