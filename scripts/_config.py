"""Load config/params.yaml with ${MTCOV_ROOT} expanded.

The scratch layout used to appear as seven absolute paths inside params.yaml,
which made a cluster-specific prefix look like a parameter of the study. It is
not: it is where the big intermediates happen to live. Writing them as
${MTCOV_ROOT}/fastq and friends leaves one variable to set, and the same one the
batch scripts already read (scripts/_env.sh).

Unset variables are left as written rather than expanded to nothing, so a
mistyped name shows up as a path containing a dollar sign instead of silently
resolving to the filesystem root.
"""

import os
import yaml

# The layout this study ran in. Kept as the default so the published
# invocations work unchanged, not because anything depends on it.
DEFAULT_ROOT = "/arf/scratch/suysal/mtcovmap"


def expand(o):
    """Expand environment variables in every string of a nested structure."""
    if isinstance(o, str):
        return os.path.expandvars(o)
    if isinstance(o, dict):
        return {k: expand(v) for k, v in o.items()}
    if isinstance(o, list):
        return [expand(v) for v in o]
    return o


def load(path):
    os.environ.setdefault("MTCOV_ROOT", DEFAULT_ROOT)
    with open(path) as f:
        return expand(yaml.safe_load(f))


def selftest():
    os.environ["MTCOV_ROOT"] = "/scratch/x"
    got = expand({"a": "${MTCOV_ROOT}/bam", "b": [1, "${MTCOV_ROOT}"], "c": 3})
    assert got == {"a": "/scratch/x/bam", "b": [1, "/scratch/x"], "c": 3}, got
    # An unset variable stays visible instead of collapsing to "" -> "/fastq".
    assert expand("${MTCOV_NO_SUCH_VAR}/fastq") == "${MTCOV_NO_SUCH_VAR}/fastq"
    # Non-strings are returned untouched, including the numeric thresholds.
    assert expand(0.05) == 0.05 and expand(None) is None
    print("ok")


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        selftest()
