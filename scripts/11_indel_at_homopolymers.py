#!/usr/bin/env python3
"""Do the poly-C tracts look noisy because indels arrive as substitutions?

The manuscript argues that at 303-315 and 16,180-16,195 an aligner may resolve a
length difference as mismatches inside the tract rather than as a gap, so a
substitution-only pileup reports them as apparent alternate alleles. That is a
claim about read alignments, and it is testable on the BAMs we already have: if
it holds, the apparent substitutions at these positions should be carried
disproportionately by reads that also carry an indel nearby.

No realignment; this reads the existing chrM BAMs.
"""
import argparse, sys
from pathlib import Path
import pandas as pd
import _config

# Tracts named in the manuscript, plus two control windows of comparable depth
# that are not homopolymeric.
TRACTS = {"poly_c_hvs2": (303, 315), "poly_c_hvs1": (16180, 16195)}
CONTROLS = {"rnr2_control": (2000, 2015), "co1_control": (6000, 6015)}
FLANK = 10          # bp either side within which an indel counts as "nearby"


def classify(bam, regions, flank=FLANK):
    """Per position: reads supporting a substitution, split by nearby-indel status."""
    import pysam
    rows = []
    with pysam.AlignmentFile(bam, "rb") as f:
        ref = next(r for r in f.references
                   if r in ("chrM", "MT", "chrM_rCRS", "NC_012920.1"))
        for name, (lo, hi) in regions.items():
            # Whether a read carries an indel is a property of the read, so it is
            # read once here and reused for every position it covers.
            indel = {}
            for read in f.fetch(ref, max(lo - flank - 1, 0), hi + flank):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                indel[read.query_name + str(read.is_read2)] = any(
                    op in (1, 2) for op, _ in (read.cigartuples or []))
            for col in f.pileup(ref, lo - 1, hi, truncate=True,
                                min_base_quality=20, stepper="nofilter"):
                pos = col.reference_pos + 1
                sub = {True: 0, False: 0}      # keyed by "read carries an indel"
                tot = {True: 0, False: 0}
                base_at = {}
                for pr in col.pileups:
                    r = pr.alignment
                    if r.is_unmapped or r.is_secondary or r.is_supplementary:
                        continue
                    key = indel.get(r.query_name + str(r.is_read2), False)
                    tot[key] += 1
                    if pr.query_position is None:        # this read has a gap here
                        continue
                    b = r.query_sequence[pr.query_position]
                    base_at[b] = base_at.get(b, 0) + 1
                    sub[key] += 0                        # filled after the major base
                # The reference base is not available here, so the majority base
                # stands in for it: these are near-fixed positions in bulk data.
                if not base_at:
                    continue
                major = max(base_at, key=base_at.get)
                for pr in col.pileups:
                    r = pr.alignment
                    if (r.is_unmapped or r.is_secondary or r.is_supplementary
                            or pr.query_position is None):
                        continue
                    if r.query_sequence[pr.query_position] != major:
                        sub[indel.get(r.query_name + str(r.is_read2), False)] += 1
                rows.append({"region": name, "pos": pos,
                             "depth_indel_reads": tot[True],
                             "depth_clean_reads": tot[False],
                             "alt_from_indel_reads": sub[True],
                             "alt_from_clean_reads": sub[False]})
    return pd.DataFrame(rows)


def summarise(d):
    """Enrichment of apparent substitutions in indel-carrying reads."""
    g = d.groupby("region")[["depth_indel_reads", "depth_clean_reads",
                             "alt_from_indel_reads", "alt_from_clean_reads"]].sum()
    g["alt_rate_indel_reads"] = g.alt_from_indel_reads / g.depth_indel_reads
    g["alt_rate_clean_reads"] = g.alt_from_clean_reads / g.depth_clean_reads
    g["enrichment"] = g.alt_rate_indel_reads / g.alt_rate_clean_reads
    # The headline quantity. The two resolutions are alternatives within one read,
    # so "substitutions enriched in indel-carrying reads" only holds where the
    # aligner gaps a minority of reads; how often it gaps any read at all is the
    # measure of length ambiguity that does not depend on which one it prefers.
    total = g.depth_indel_reads + g.depth_clean_reads
    g["frac_reads_with_indel"] = g.depth_indel_reads / total
    g["overall_alt_rate"] = (g.alt_from_indel_reads + g.alt_from_clean_reads) / total
    g["frac_alt_carried_by_indel_reads"] = (
        g.alt_from_indel_reads / (g.alt_from_indel_reads + g.alt_from_clean_reads))
    return g.reset_index()


def selftest():
    d = pd.DataFrame([{"region": "t", "pos": 1, "depth_indel_reads": 100,
                       "depth_clean_reads": 100, "alt_from_indel_reads": 50,
                       "alt_from_clean_reads": 5}])
    s = summarise(d).iloc[0]
    assert abs(s.enrichment - 10.0) < 1e-9
    assert abs(s.frac_reads_with_indel - 0.5) < 1e-9
    assert abs(s.frac_alt_carried_by_indel_reads - 50 / 55) < 1e-9
    print("selftest ok")


def main():
    a = argparse.ArgumentParser()
    a.add_argument("--config", default="config/params.yaml")
    a.add_argument("--samples", default="config/samples.tsv")
    a.add_argument("--out", default="results/tables")
    a.add_argument("--arm", default="A")
    a.add_argument("--selftest", action="store_true")
    a = a.parse_args()
    if a.selftest:
        return selftest()
    c = _config.load(a.config)
    bam_dir = Path(c["alignment"]["bam_dir"])
    samples = pd.read_csv(a.samples, sep="\t")
    regions = {**TRACTS, **CONTROLS}
    per = []
    for r in samples.itertuples():
        bam = bam_dir / f"{r.run_accession}_{a.arm}_chrM.bam"
        if not bam.exists():
            print(f"missing {bam}", file=sys.stderr); continue
        d = classify(bam, regions)
        per.append(d.assign(run_accession=r.run_accession, tissue=r.tissue))
        print(f"  {r.run_accession}", file=sys.stderr)
    d = pd.concat(per, ignore_index=True)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    d.to_csv(out / "homopolymer_indel_per_position.tsv", sep="\t", index=False)
    s = summarise(d)
    s.to_csv(out / "homopolymer_indel_summary.tsv", sep="\t", index=False)
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
