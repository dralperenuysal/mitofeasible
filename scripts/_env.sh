# Shared environment for every batch script. Sourced, never executed.
#
# The point of this file is that the same seven lines used to be copied into
# each .slurm, which made the cluster assumptions look structural when they are
# not. Everything below has an override, so the same scripts run under sbatch on
# TRUBA and in a plain loop on a workstation:
#
#   MTCOV_ROOT   scratch: FASTQ, BAMs, coverage, the container   (big, not in
#                git) -- required by every phase that touches data, and with no
#                default, because there is no sensible one
#   MTCOV_REPO   this checkout -- required, and the one thing with no default:
#                a batch script cannot find itself, because SLURM copies it to
#                the node's spool before running it, so $BASH_SOURCE points at
#                the copy. run_all.sh sets it from its own location and exports
#                it, and sbatch passes it on; only hand-submission must export
#                it first, and then it says so rather than guessing.
#   MTCOV_SIF    the Apptainer image
#   MTCOV_BIND   paths to bind into the container
#
# Off the cluster, set MTCOV_ROOT and MTCOV_REPO; the rest follows.

ROOT=${MTCOV_ROOT:-}                               # scratch: big intermediates
REPO=${MTCOV_REPO:?set MTCOV_REPO to this checkout}  # this checkout
SIF=${MTCOV_SIF:-${ROOT:+$ROOT/mtcovmap.sif}}
# Exported so the Python side sees the same layout: config/params.yaml writes
# every directory as ${MTCOV_ROOT}/... and scripts/_config.py expands it.
# ROOT has no default -- there is no sensible one, and the previous default was
# one person's scratch directory, which every other user would have silently
# inherited. Unset, the phases that need it say so; plan and selftest do not.
export MTCOV_REPO="$REPO"
[ -n "$ROOT" ] && export MTCOV_ROOT="$ROOT"

# Job logs are written to logs/ relative to the submit directory, not to an
# absolute path. SLURM opens those files before the script body runs, so an
# absolute path that does not exist on the cluster kills the job before any
# mkdir can help -- and leaves nothing in any log to explain it. Relative keeps
# the same destination here (run_all.sh chdirs to ROOT first) and works
# unchanged elsewhere. Whoever submits by hand must chdir to ROOT too.
#
# Best-effort: `run_all.sh plan` and `selftest` are meant to work in a fresh
# clone with no scratch anywhere, so an unreachable ROOT is not an error here.
# The phases that need it call need_root, which is where it is an error.
[ -n "$ROOT" ] && { mkdir -p "$ROOT/logs" 2>/dev/null || true; }

# Bind /arf when it exists, because on TRUBA both ROOT and REPO live under it
# and one bind covers both. Elsewhere bind exactly the two directories used.
BIND=${MTCOV_BIND:-$([ -d /arf ] && echo /arf || echo "${ROOT:-$REPO},$REPO")}

# The container is how the published numbers were produced, so it is the default
# wherever it can be used. Where apptainer is absent the scripts still run
# against whatever interpreter is on PATH -- with the versions unpinned, which
# is the user's problem to solve and not a reason to refuse to start.
if command -v apptainer >/dev/null 2>&1 && [ -n "$SIF" ] && [ -f "$SIF" ]; then
    APPT="apptainer exec --bind $BIND $SIF"
else
    APPT=""
    echo "note: no container${SIF:+ at $SIF}; running on PATH (versions unpinned)" >&2
fi

# Word-splitting is deliberate: APPT is a command plus flags, and quoting it
# would make the whole string one command name. Use `pyrun` rather than $APPT
# directly so the empty case collapses to a bare interpreter.
pyrun() { if [ -n "$APPT" ]; then $APPT python3 "$@"; else python3 "$@"; fi; }
run_in() { if [ -n "$APPT" ]; then $APPT "$@"; else "$@"; fi; }
