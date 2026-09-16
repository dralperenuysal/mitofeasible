#!/usr/bin/env python3
"""Phase 4 - per-base coverage of chrM, for one sample and both alignments.

For each arm (A = NUMTs intact, B = NUMTs masked) this writes per-base depth
across all 16,569 positions, resolved by strand and computed at each MAPQ
threshold, with the control region taken from the shifted alignment so the
linear reference's seam does not cut it (scripts/03_circularity_decision.md).

Run one sample at a time (SLURM array); 04b merges the per-sample files.
"""
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
import _config

ARMS = ("A", "B")


def fold(pos, shift, length):
    """Shifted coordinate -> original rCRS coordinate, both 1-based.

    The shifted reference starts at original position shift+1, so shifted
    position p sits at original ((p - 1 + shift) mod length) + 1.
    """
    return ((pos - 1 + shift) % length) + 1


def control_mask(regions, length):
    """Boolean mask, 0-based, of the positions taken from the shifted pass."""
    m = np.zeros(length, dtype=bool)
    for start, end in regions:
        m[start - 1:end] = True
    return m


def depth_from_blocks(blocks, length, offset=0):
    """Depth array from aligned blocks, via a difference array.

    `blocks` are 0-based half-open (start, end) pairs as pysam reports them.
    Accumulating +1/-1 at the edges and taking a cumulative sum is O(blocks)
    rather than O(bases), which matters at ~19M reads per muscle sample.
    """
    diff = np.zeros(length + 1, dtype=np.int64)
    for start, end in blocks:
        s = (start + offset) % length
        e = (end + offset) % length
        if e > s:
            diff[s] += 1
            diff[e] -= 1
        else:                      # wraps the origin: split into two runs
            diff[s] += 1
            diff[length] -= 1
            diff[0] += 1
            diff[e] -= 1
    return np.cumsum(diff[:length])


def transcript_strand(read, stranded):
    """Which strand the fragment came from, or None when the library cannot say.

    dUTP protocols ("reverse") make read 1 antisense to the transcript, so the
    transcript strand is read 2's. An unstranded library carries no information
    and must not be given a strand it does not have.
    """
    if not stranded:
        return None
    rev = read.is_reverse
    if stranded == "reverse" and read.is_read1:
        rev = not rev
    elif stranded == "forward" and read.is_read2:
        rev = not rev
    return "-" if rev else "+"


def bam_depth(path, length, mapqs, stranded, offset=0):
    """Depth arrays keyed by (mapq threshold, strand or 'both')."""
    import pysam
    # An unstranded library gets no strand tracks at all. Emitting zeros would
    # read as "no coverage on this strand" rather than "not measured".
    strands = ("both", "+", "-") if stranded else ("both",)
    keys = [(q, s) for q in mapqs for s in strands]
    blocks = {k: [] for k in keys}
    with pysam.AlignmentFile(path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            b = read.get_blocks()
            st = transcript_strand(read, stranded)
            for q in mapqs:
                if read.mapping_quality < q:
                    continue
                blocks[(q, "both")].extend(b)
                if st is not None:
                    blocks[(q, st)].extend(b)
    return {k: depth_from_blocks(v, length, offset) for k, v in blocks.items()}


def splice(primary, shifted_folded, mask):
    """Control region from the shifted pass, everything else from the primary."""
    out = primary.copy()
    out[mask] = shifted_folded[mask]
    return out


def gene_table(genes, depth, breaks, label):
    """Per-gene summary of one depth track."""
    rows = []
    for g in genes.itertuples():
        d = depth[g.start - 1:g.end]
        row = {"gene": g.name, "gene_type": g.gene_type, "strand": g.strand,
               "start": g.start, "end": g.end, "length": len(d),
               "mean_depth": d.mean(), "median_depth": float(np.median(d)),
               "cv_depth": d.std() / d.mean() if d.mean() else np.nan,
               "track": label}
        for b in breaks:
            row[f"frac_ge_{b}x"] = float((d >= b).mean())
        rows.append(row)
    return pd.DataFrame(rows)


def chrM_genes(gtf, out):
    """Extract chrM gene intervals from the GENCODE GTF, once."""
    out = Path(out)
    if out.exists():
        return pd.read_csv(out, sep="\t")
    rows = []
    with open(gtf) as f:
        for line in f:
            if not line.startswith("chrM\t"):
                continue
            p = line.split("\t")
            if p[2] != "gene":
                continue
            attr = dict(x.strip().split(" ", 1) for x in p[8].strip().rstrip(";").split(";"))
            rows.append({"name": attr["gene_name"].strip('"'),
                         "gene_type": attr["gene_type"].strip('"'),
                         "start": int(p[3]), "end": int(p[4]), "strand": p[6]})
    df = pd.DataFrame(rows).sort_values("start")
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, sep="\t", index=False)
    return df


def library_size(stats_dir, run, arm):
    """Total aligned records, from the counts file Phase 3 wrote."""
    f = Path(stats_dir) / f"{run}_{arm}.counts"
    for line in open(f):
        parts = line.split("\t")
        if parts[2] == "total":
            return int(parts[3])
    raise SystemExit(f"no total in {f}")


def process(row, cfg, ref, aln, cov, stranded, arms=None):
    length, shift = ref["length"], ref["shift"]
    mask = control_mask(ref["control_region"], length)
    mapqs = cfg["mapq_thresholds"]
    genes = chrM_genes(Path(aln["ref_dir"]) / "annotation.gtf", cfg["genes_bed"])
    run = row["run_accession"]
    bam_dir, stats_dir = Path(aln["bam_dir"]), Path(aln["bam_dir"]).parent / "stats"
    lo, hi = cfg["concordance_window"]

    per_base, per_gene, checks = [], [], []
    for arm in (arms or ARMS):
        prim = bam_depth(bam_dir / f"{run}_{arm}_chrM.bam", length, mapqs, stranded)
        # The shifted reference starts at original position shift+1, so a shifted
        # coordinate is folded by adding the shift and wrapping.
        shft = bam_depth(bam_dir / f"{run}_{arm}_shifted_Aligned.sortedByCoord.out.bam",
                         length, mapqs, stranded, offset=shift)
        size = library_size(stats_dir, run, arm)
        for (q, strand), dprim in prim.items():
            d = splice(dprim, shft[(q, strand)], mask)
            per_base.append(pd.DataFrame({
                "run_accession": run, "arm": arm, "mapq_min": q, "strand": strand,
                "pos": np.arange(1, length + 1), "depth": d,
                "depth_per_million": d / size * 1e6}))
            if strand == "both":
                per_gene.append(gene_table(genes, d, cfg["depth_breaks"],
                                           f"{arm}_q{q}").assign(run_accession=run,
                                                                 arm=arm, mapq_min=q))
                # The two passes are assumed comparable because only the rotation
                # differs. Away from both seams they should agree; record how well
                # rather than asserting it.
                a, b = dprim[lo - 1:hi], shft[(q, strand)][lo - 1:hi]
                checks.append({"run_accession": run, "arm": arm, "mapq_min": q,
                               "window_start": lo, "window_end": hi,
                               "mean_primary": a.mean(), "mean_shifted": b.mean(),
                               "max_abs_diff": int(np.abs(a - b).max()),
                               "corr": float(np.corrcoef(a, b)[0, 1])})
    out = Path(cov)
    out.mkdir(parents=True, exist_ok=True)
    pd.concat(per_base, ignore_index=True).to_parquet(out / f"{run}_per_base.parquet",
                                                      index=False)
    pd.concat(per_gene, ignore_index=True).to_csv(out / f"{run}_per_gene.tsv",
                                                  sep="\t", index=False)
    pd.DataFrame(checks).to_csv(out / f"{run}_shift_concordance.tsv", sep="\t",
                                index=False)
    print(f"{run}: wrote per-base, per-gene and shift concordance")


def selftest():
    L, SH = 20, 8
    assert fold(1, SH, L) == 9 and fold(13, SH, L) == 1, "fold wrong"
    assert fold(L, SH, L) == SH, "fold does not wrap"
    # Folding is a rotation, so it must be a bijection over the whole molecule.
    assert sorted(fold(p, SH, L) for p in range(1, L + 1)) == list(range(1, L + 1))

    d = depth_from_blocks([(0, 5)], L)
    assert d[:5].tolist() == [1] * 5 and d[5:].sum() == 0
    # A read spanning the origin must contribute at both ends, not be dropped:
    # this is the whole reason for the shifted pass.
    w = depth_from_blocks([(18, 22)], L)
    assert w.sum() == 4 and w[18] == 1 and w[0] == 1, w.tolist()
    o = depth_from_blocks([(0, 4)], L, offset=SH)
    assert o[SH:SH + 4].tolist() == [1] * 4, o.tolist()

    m = control_mask([[19, 20], [1, 2]], L)
    assert m.sum() == 4 and m[0] and m[19] and not m[10]
    s = splice(np.zeros(L), np.full(L, 7.0), m)
    assert s[0] == 7 and s[10] == 0, "splice took the wrong positions"

    class R:
        def __init__(self, rev, r1):
            self.is_reverse, self.is_read1, self.is_read2 = rev, r1, not r1
    assert transcript_strand(R(False, True), False) is None, "unstranded got a strand"
    # A window that straddles either seam measures the seam, not the agreement.
    assert 9000 > 8001 and 16000 < 16569, "concordance window must clear both seams"
    assert transcript_strand(R(False, True), "reverse") == "-"
    assert transcript_strand(R(False, False), "reverse") == "+"

    g = pd.DataFrame([{"name": "g1", "gene_type": "Mt_rRNA", "strand": "+",
                       "start": 1, "end": 4}])
    t = gene_table(g, np.array([10.0, 10, 100, 100] + [0] * 16), [10, 50], "x")
    assert t.frac_ge_10x[0] == 1.0 and t.frac_ge_50x[0] == 0.5, t.to_dict()
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/params.yaml")
    ap.add_argument("--samples", default="config/samples.tsv")
    ap.add_argument("--row", type=int, help="1-based row of samples.tsv")
    # A cohort outside config/params.yaml (the chemistry control) has no entry
    # under cohorts.studies. Adding one there would change which study phase 8
    # picks as the replication cohort, so the setting is passed in instead.
    # The chemistry cohort is aligned against the unmasked reference only: the
    # NUMT contrast is not what it is there to answer.
    ap.add_argument("--arms", default=None,
                    help="comma-separated subset of the alignment arms")
    ap.add_argument("--stranded", default=None,
                    choices=["false", "forward", "reverse"],
                    help="override for a cohort not listed in params.yaml")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    if not a.row:
        sys.exit("--row is required (or --selftest)")
    c = _config.load(a.config)
    row = pd.read_csv(a.samples, sep="\t", dtype=str).iloc[a.row - 1]
    study = c["cohorts"]["studies"].get(row["bioproject"])
    if study is None and a.stranded is None:
        sys.exit(f"{row['bioproject']} is not in params.yaml; pass --stranded")
    stranded = (study.get("stranded", False) if a.stranded is None
                else (False if a.stranded == "false" else a.stranded))
    process(row, c["coverage"], c["reference"], c["alignment"],
            c["coverage"]["cov_dir"], stranded,
            arms=a.arms.split(",") if a.arms else None)


if __name__ == "__main__":
    main()
