#!/usr/bin/env python3
"""MitoFeasible - what can your own RNA-seq actually detect on chrM?

Point it at a chrM BAM and it returns, per position, the lowest allele fraction
a variant call there is supported at:

    from mitofeasible import feasibility
    m = feasibility("sample.bam")
    m.loc[4500]
    # depth 47, min_detectable_af 0.112  -> a 3% signal there is noise

or from the shell:

    python mitofeasible.py sample.bam -o feasibility.tsv
    python mitofeasible.py sample.bam --calls my_calls.tsv

The transferable half of the model is the per-position error rate measured in
this study, which a user cannot get from their own data: at a position where the
donor's DNA carries no variant, the non-reference reads in their RNA are error -
sequencing error, reverse-transcription misincorporation at modified RNA bases,
and NUMT-derived reads. The other half, depth, is read from the user's BAM,
because depth does not transfer between datasets at all.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

CHRM_LENGTH = 16569
DEFAULT_ERRORS = Path(__file__).resolve().parent.parent / "results" / "tables" / "error_rate_per_position.tsv"
DEFAULT_ALPHA = 0.01
DEFAULT_N_TESTS = 16569
DEFAULT_POWER = 0.8


def min_alt_reads(depth, eps, alpha):
    """Smallest alt count not explainable by error at significance `alpha`.

    One-sided: the question is only whether the alt load is too high for
    Binomial(depth, eps).
    """
    # Floor, not round: depth is a count, and a median over an even number of
    # samples is not one. scipy's binom returns NaN for non-integer n without
    # raising, which silently turns every downstream comparison False.
    d = np.floor(np.atleast_1d(np.asarray(depth, dtype=float)))
    k = stats.binom.isf(alpha, d, eps) + 1
    return np.maximum(k, 1)


def lod(depth, eps, alpha, power, tol=1e-6):
    """Limit of detection: smallest true allele fraction found with `power`.

    A threshold alone says whether a count is surprising. It does not say whether
    a variant at frequency f would have produced such a count often enough to be
    found. Only the second is useful to someone planning an analysis.
    """
    d = np.floor(np.atleast_1d(np.asarray(depth, dtype=float)))
    k = min_alt_reads(d, eps, alpha)
    lo, hi = np.zeros_like(d), np.ones_like(d)
    for _ in range(60):
        mid = (lo + hi) / 2
        ok = stats.binom.sf(k - 1, d, mid) >= power
        hi = np.where(ok, mid, hi)
        lo = np.where(ok, lo, mid)
        if np.all(hi - lo < tol):
            break
    out = np.where(d > 0, hi, np.nan)
    return out if out.size > 1 else float(out[0])


def chrm_depth(bam, min_base_quality=20, min_mapping_quality=0, length=CHRM_LENGTH):
    """Per-base depth on chrM, 1-based, from any BAM that contains it."""
    import pysam
    depth = np.zeros(length, dtype=np.int64)
    with pysam.AlignmentFile(bam, "rb") as b:
        names = set(b.references)
        # Callers name the same sequence chrM, MT or NC_012920.1 depending on the
        # reference they used; refusing any of them would be a pointless failure.
        ref = next((n for n in ("chrM", "MT", "chrM_rCRS", "NC_012920.1") if n in names), None)
        if ref is None:
            raise ValueError(f"no chrM-like sequence in {bam}; found {sorted(names)[:8]}")
        if b.get_reference_length(ref) != length:
            raise ValueError(f"{ref} is {b.get_reference_length(ref)} bp, expected {length};"
                             " the map is defined on rCRS (NC_012920.1)")
        for col in b.pileup(ref, max_depth=1_000_000, min_base_quality=min_base_quality,
                            min_mapping_quality=min_mapping_quality, stepper="samtools"):
            if col.reference_pos < length:
                depth[col.reference_pos] = col.get_num_aligned()
    return depth


def load_errors(path=None):
    e = pd.read_csv(path or DEFAULT_ERRORS, sep="\t")
    return e[["pos", "epsilon"]]


def feasibility(bam, errors=None, alpha=DEFAULT_ALPHA, n_tests=DEFAULT_N_TESTS,
                power=DEFAULT_POWER, **pileup_kw):
    """Per-position feasibility for one BAM. Returns a frame indexed by position."""
    depth = chrm_depth(bam, **pileup_kw)
    e = load_errors(errors)
    df = pd.DataFrame({"pos": np.arange(1, CHRM_LENGTH + 1), "depth": depth})
    df = df.merge(e, on="pos", how="left")
    # Positions the study could not measure fall back to its median rate, and say
    # so, rather than dropping out of the map without explanation.
    df["epsilon_imputed"] = df.epsilon.isna()
    df["epsilon"] = df.epsilon.fillna(e.epsilon.median())
    df["min_detectable_af"] = lod(df.depth.values, df.epsilon.values,
                                  alpha / n_tests, power)
    return df.set_index("pos")


def check_calls(calls, fmap):
    """Flag calls that sit below what their position supports.

    `calls` needs a `pos` and an allele-fraction column (`af` or `alt_frac`).
    """
    c = calls.copy()
    af = "af" if "af" in c.columns else "alt_frac"
    j = c.join(fmap[["depth", "epsilon", "min_detectable_af"]], on="pos")
    j["supported"] = j[af] >= j.min_detectable_af
    return j


def _selftest():
    assert lod(10000, 1e-3, 1e-6, .8) < lod(100, 1e-3, 1e-6, .8), "deeper must be easier"
    assert lod(1000, 1e-2, 1e-6, .8) > lod(1000, 1e-4, 1e-6, .8), "noisier must be harder"
    assert np.isnan(lod(0, 1e-3, 1e-6, .8)), "zero depth is not detectable"
    half = lod(477.5, 4.68e-4, 1e-6, .8)
    assert np.isfinite(half) and abs(half - lod(477, 4.68e-4, 1e-6, .8)) < 1e-12

    fm = pd.DataFrame({"pos": [10, 20], "depth": [50, 20000], "epsilon": [1e-3, 1e-3],
                       "min_detectable_af": [0.20, 0.002]}).set_index("pos")
    calls = pd.DataFrame({"pos": [10, 20], "af": [0.05, 0.05]})
    out = check_calls(calls, fm)
    assert out.supported.tolist() == [False, True], out.to_dict("list")
    print("selftest ok")


def main():
    import argparse
    ap = argparse.ArgumentParser(description="Per-position detection limits for a chrM BAM")
    ap.add_argument("bam", nargs="?", help="indexed BAM containing chrM")
    ap.add_argument("-o", "--out", help="write the map here (TSV); default stdout summary")
    ap.add_argument("--calls", help="TSV with pos and af columns; flags unsupported calls")
    ap.add_argument("--errors", help="error-rate table (default: the one shipped with this study)")
    ap.add_argument("--alpha", type=float, default=DEFAULT_ALPHA)
    ap.add_argument("--power", type=float, default=DEFAULT_POWER)
    ap.add_argument("--min-mapping-quality", type=int, default=0)
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return _selftest()
    if not a.bam:
        ap.error("a BAM is required (or --selftest)")
    m = feasibility(a.bam, errors=a.errors, alpha=a.alpha, power=a.power,
                    min_mapping_quality=a.min_mapping_quality)
    if a.out:
        m.to_csv(a.out, sep="\t")
        print(f"wrote {a.out}")
    if a.calls:
        calls = pd.read_csv(a.calls, sep="\t")
        j = check_calls(calls, m)
        print(j.to_string(index=False))
        print(f"\n{int((~j.supported).sum())} of {len(j)} calls fall below the "
              "detection limit at their position")
    if not a.out and not a.calls:
        q = m.min_detectable_af.quantile([.5, .9, .99])
        print(f"median depth {int(m.depth.median())}x")
        print(f"detectable AF: median {q[.5]*100:.2f}%  p90 {q[.9]*100:.2f}%  p99 {q[.99]*100:.2f}%")
        for t in (0.01, 0.02, 0.05):
            print(f"positions supporting {int(t*100)}%: "
                  f"{(m.min_detectable_af <= t).mean()*100:.1f}%")


if __name__ == "__main__":
    main()
