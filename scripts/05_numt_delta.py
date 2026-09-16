#!/usr/bin/env python3
"""Phase 5 - per-position effect of NUMT masking. Primary result.

Alignment A leaves the nuclear NUMTs in place, B masks them. Every other
parameter is identical, so the per-position difference B - A is the depth that
NUMTs were taking from chrM, position by position.

Cohorts are summarised separately and never pooled: they differ in tissue and in
study, and the replication cohort is not evidence for the primary one until
Phase 8 (AGENTS.md 7.3).
"""
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
import _config


def delta_table(per_base):
    """Wide per-position A/B depths -> signed and relative deltas.

    Normalised depth is used, not raw: A and B have slightly different library
    sizes because masking removes alignments elsewhere in the genome, and the
    raw difference would carry that along.
    """
    idx = ["run_accession", "mapq_min", "strand", "pos"]
    w = per_base.pivot_table(index=idx, columns="arm",
                             values=["depth", "depth_per_million"]).reset_index()
    w.columns = [c[0] if not c[1] else f"{c[0]}_{c[1]}" for c in w.columns]
    w["delta"] = w["depth_per_million_B"] - w["depth_per_million_A"]
    # Relative to A, the alignment a reader would produce without masking. Where
    # A is zero the ratio is undefined rather than infinite; such positions are
    # reported through the absolute delta instead.
    a = w["depth_per_million_A"]
    w["rel_delta"] = np.where(a > 0, w["delta"] / a.where(a > 0), np.nan)
    return w


def summarise_positions(deltas, samples, thresholds, min_frac):
    """Per (cohort, mapq, strand, position) summary across that cohort's samples."""
    d = deltas.merge(samples[["run_accession", "tissue", "role"]], on="run_accession")
    keys = ["tissue", "role", "mapq_min", "strand", "pos"]
    g = d.groupby(keys)
    out = g.agg(n_samples=("run_accession", "nunique"),
                mean_depth_A=("depth_per_million_A", "mean"),
                mean_depth_B=("depth_per_million_B", "mean"),
                mean_delta=("delta", "mean"),
                median_rel_delta=("rel_delta", "median"),
                max_rel_delta=("rel_delta", "max")).reset_index()
    for t in thresholds:
        frac = g["rel_delta"].apply(lambda s, t=t: (s.abs() >= t).mean())
        out[f"frac_samples_ge_{t}"] = frac.values
        out[f"affected_{t}"] = out[f"frac_samples_ge_{t}"] >= min_frac
    return out


def call_regions(flags, positions, min_length):
    """Contiguous runs of affected positions, as (start, end) 1-based inclusive.

    Runs shorter than min_length are dropped: one position can move on a single
    read's placement, a run cannot.
    """
    flags = np.asarray(flags, dtype=bool)
    positions = np.asarray(positions)
    order = np.argsort(positions)
    flags, positions = flags[order], positions[order]
    regions, start = [], None
    for i, f in enumerate(flags):
        if f and start is None:
            start = positions[i]
        elif not f and start is not None:
            if positions[i - 1] - start + 1 >= min_length:
                regions.append((int(start), int(positions[i - 1])))
            start = None
    if start is not None and positions[-1] - start + 1 >= min_length:
        regions.append((int(start), int(positions[-1])))
    return regions


def annotate(regions, genes):
    """Name each region by the genes it overlaps, or 'intergenic'."""
    rows = []
    for start, end in regions:
        hit = genes[(genes.start <= end) & (genes.end >= start)]
        rows.append({"start": start, "end": end, "length": end - start + 1,
                     "genes": ",".join(hit.name) if len(hit) else "intergenic",
                     "gene_types": ",".join(sorted(set(hit.gene_type))) if len(hit)
                     else "intergenic"})
    return pd.DataFrame(rows)


def run(cov_dir, samples, genes, cfg, out_dir):
    files = sorted(Path(cov_dir).glob("*_per_base.parquet"))
    if not files:
        sys.exit(f"no per-base coverage in {cov_dir}")
    print(f"reading {len(files)} coverage files", file=sys.stderr)
    deltas = pd.concat([delta_table(pd.read_parquet(f)) for f in files],
                       ignore_index=True)
    pos = summarise_positions(deltas, samples, cfg["rel_thresholds"],
                              cfg["min_sample_fraction"])

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    pos.to_parquet(out_dir / "numt_delta_per_position.parquet", index=False)

    region_rows = []
    for t in cfg["rel_thresholds"]:
        for (tissue, role, q, strand), sub in pos.groupby(["tissue", "role",
                                                           "mapq_min", "strand"]):
            regs = call_regions(sub[f"affected_{t}"], sub["pos"],
                                cfg["min_region_length"])
            if not regs:
                continue
            df = annotate(regs, genes)
            df.insert(0, "strand", strand)
            df.insert(0, "mapq_min", q)
            df.insert(0, "role", role)
            df.insert(0, "tissue", tissue)
            df.insert(0, "rel_threshold", t)
            region_rows.append(df)
    regions = (pd.concat(region_rows, ignore_index=True) if region_rows
               else pd.DataFrame(columns=["rel_threshold", "tissue", "start", "end"]))
    regions.to_csv(out_dir / "numt_delta_regions.tsv", sep="\t", index=False)

    # Per-sample totals: the cohort-level number a reader wants first.
    per_sample = (deltas[deltas.strand == "both"]
                  .groupby(["run_accession", "mapq_min"])
                  .agg(mean_rel_delta=("rel_delta", "mean"),
                       max_rel_delta=("rel_delta", "max"),
                       positions_ge_5pct=("rel_delta",
                                          lambda s: int((s.abs() >= 0.05).sum())))
                  .reset_index()
                  .merge(samples[["run_accession", "tissue", "role"]],
                         on="run_accession"))
    per_sample.to_csv(out_dir / "numt_delta_per_sample.tsv", sep="\t", index=False)
    print(f"wrote per-position, region and per-sample tables to {out_dir}")
    return pos, regions, per_sample


def selftest():
    pb = pd.DataFrame({
        "run_accession": ["r1"] * 6, "arm": ["A", "B"] * 3, "mapq_min": [0] * 6,
        "strand": ["both"] * 6, "pos": [1, 1, 2, 2, 3, 3],
        "depth": [100, 110, 50, 50, 0, 5],
        "depth_per_million": [10.0, 11.0, 5.0, 5.0, 0.0, 0.5]})
    w = delta_table(pb).sort_values("pos")
    assert list(w.delta) == [1.0, 0.0, 0.5], list(w.delta)
    assert abs(w.rel_delta.iloc[0] - 0.1) < 1e-9
    # A zero-depth baseline must not become an infinite ratio.
    assert np.isnan(w.rel_delta.iloc[2]), w.rel_delta.tolist()

    # A run shorter than min_length is not a region; the long one is.
    flags = [False, True, False] + [True] * 5 + [False]
    pos = list(range(1, 10))
    assert call_regions(flags, pos, 5) == [(4, 8)], call_regions(flags, pos, 5)
    assert call_regions(flags, pos, 1) == [(2, 2), (4, 8)]
    # A run touching the end of the array must still close.
    assert call_regions([False, True, True], [1, 2, 3], 2) == [(2, 3)]

    genes = pd.DataFrame([{"name": "MT-CO1", "gene_type": "protein_coding",
                           "start": 5, "end": 10},
                          {"name": "MT-TS1", "gene_type": "Mt_tRNA",
                           "start": 50, "end": 60}])
    a = annotate([(6, 8), (20, 30)], genes)
    assert a.genes.tolist() == ["MT-CO1", "intergenic"], a.genes.tolist()

    d = delta_table(pb).assign(rel_delta=[0.10, 0.0, np.nan])
    s = pd.DataFrame({"run_accession": ["r1"], "tissue": ["lcl"], "role": ["primary"]})
    out = summarise_positions(d, s, [0.05], 0.5).sort_values("pos")
    assert out["affected_0.05"].tolist() == [True, False, False], out.to_dict("list")
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/params.yaml")
    ap.add_argument("--samples", default="config/samples.tsv")
    ap.add_argument("--out", default="results/tables")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    c = _config.load(a.config)
    samples = pd.read_csv(a.samples, sep="\t")
    genes = pd.read_csv(c["coverage"]["genes_bed"], sep="\t")
    run(c["coverage"]["cov_dir"], samples, genes, c["numt_delta"], a.out)


if __name__ == "__main__":
    main()
