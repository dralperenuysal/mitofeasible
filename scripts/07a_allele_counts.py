#!/usr/bin/env python3
"""Phase 7a - per-position allele counts, the input the detection model needs.

Phase 4 measured how deep each position is. The feasibility map also needs to
know what the reads there actually say: at a position where DNA carries no
variant, whatever fraction of reads disagrees with the reference is the error
floor - sequencing error, reverse-transcription misincorporation at modified
RNA bases, and NUMT-derived reads together.

One sample per task; 07 estimates the error rate from these and builds the map.
"""
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

BASES = "ACGT"


def read_reference(path):
    return "".join(l.strip() for l in open(path) if not l.startswith(">")).upper()


def count_alleles(bam, ref, min_bq, min_mapq=0):
    """Per-position A/C/G/T counts as an (L, 4) array of int32."""
    import pysam
    L = len(ref)
    arr = np.zeros((L, 4), dtype=np.int32)
    idx = {b: i for i, b in enumerate(BASES)}
    with pysam.AlignmentFile(bam, "rb") as b:
        for col in b.pileup(max_depth=1_000_000, min_base_quality=min_bq,
                            min_mapping_quality=min_mapq, stepper="samtools"):
            pos = col.reference_pos
            if pos >= L:
                continue
            for base in col.get_query_sequences(add_indels=False):
                base = base.upper()
                j = idx.get(base)
                if j is not None:
                    arr[pos, j] += 1
    return arr


def to_frame(arr, ref, run, arm, skip):
    """Counts -> long table with the reference base, depth and non-reference load."""
    L = len(ref)
    depth = arr.sum(axis=1)
    ref_idx = np.array([BASES.find(b) for b in ref])
    # rCRS carries characters no read can match (the N placeholder at 3107). Those
    # positions get no reference count rather than a silently wrong one.
    valid = ref_idx >= 0
    ref_count = np.zeros(L, dtype=np.int64)
    ref_count[valid] = arr[np.arange(L)[valid], ref_idx[valid]]
    alt = depth - ref_count
    # The largest non-reference allele, which is what a caller would report.
    a2 = arr.copy()
    a2[np.arange(L)[valid], ref_idx[valid]] = -1
    top_alt_idx = a2.argmax(axis=1)
    top_alt = a2[np.arange(L), top_alt_idx]
    top_alt = np.where(top_alt < 0, 0, top_alt)
    df = pd.DataFrame({
        "run_accession": run, "arm": arm, "pos": np.arange(1, L + 1),
        "ref_base": list(ref), "depth": depth,
        "A": arr[:, 0], "C": arr[:, 1], "G": arr[:, 2], "T": arr[:, 3],
        "ref_count": ref_count, "alt_count": alt,
        "top_alt_base": [BASES[i] for i in top_alt_idx],
        "top_alt_count": top_alt,
    })
    df["alt_frac"] = np.where(df.depth > 0, df.alt_count / df.depth.where(df.depth > 0), np.nan)
    df["top_alt_frac"] = np.where(df.depth > 0, df.top_alt_count / df.depth.where(df.depth > 0), np.nan)
    df["usable"] = valid & ~df.pos.isin(skip).values
    return df


def selftest():
    ref = "ACGN"
    arr = np.array([[90, 10, 0, 0],     # ref A, 10% C  -> alt_frac .10
                    [0, 0, 0, 0],       # no coverage
                    [1, 2, 97, 0],      # ref G, 3% split between A and C
                    [0, 0, 0, 50]], dtype=np.int32)   # ref N: nothing can match
    d = to_frame(arr, ref, "r1", "A", [2])
    assert d.depth.tolist() == [100, 0, 100, 50]
    assert abs(d.alt_frac[0] - 0.10) < 1e-9, d.alt_frac.tolist()
    assert np.isnan(d.alt_frac[1]), "zero depth must be nan, not zero"
    assert d.top_alt_base[0] == "C" and d.top_alt_count[0] == 10
    # At the N placeholder every read counts as non-reference; the position is
    # marked unusable rather than reported as a 100% variant.
    assert d.alt_count[3] == 50 and not d.usable[3]
    assert not d.usable[1] is None
    assert d.usable.tolist() == [True, False, True, False], d.usable.tolist()
    assert d.top_alt_base[2] == "C" and d.top_alt_count[2] == 2, "largest non-ref allele"
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/params.yaml")
    ap.add_argument("--samples", default="config/samples.tsv")
    ap.add_argument("--row", type=int)
    # The chemistry cohort is aligned against the unmasked reference only.
    ap.add_argument("--arms", default="A,B",
                    help="comma-separated alignment arms to count")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if not a.row:
        sys.exit("--row is required")
    c = yaml.safe_load(open(a.config))
    acfg, hcfg = c["alignment"], c["haplogroup"]
    out = Path(c["feasibility"]["counts_dir"])
    out.mkdir(parents=True, exist_ok=True)
    ref = read_reference(Path(acfg["ref_dir"]) / "chrM_rCRS.fa")
    row = pd.read_csv(a.samples, sep="\t", dtype=str).iloc[a.row - 1]
    run = row["run_accession"]
    frames = []
    for arm in a.arms.split(","):
        arr = count_alleles(Path(acfg["bam_dir"]) / f"{run}_{arm}_chrM.bam", ref,
                            hcfg["min_base_quality"])
        frames.append(to_frame(arr, ref, run, arm, hcfg["skip_positions"]))
        print(f"{run} {arm}: median depth {int(frames[-1].depth.median())}", file=sys.stderr)
    pd.concat(frames, ignore_index=True).to_parquet(out / f"{run}_alleles.parquet",
                                                    index=False)


if __name__ == "__main__":
    main()
