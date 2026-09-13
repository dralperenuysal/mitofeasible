#!/usr/bin/env python3
"""Phase 7 - the feasibility map. Core deliverable.

For every position, the lowest allele fraction a caller could actually detect,
given how deep that position is and how noisy it is.

The error rate is measured, not assumed. At a position where a donor's own 1000
Genomes DNA carries no variant, every non-reference read in that donor's RNA is
error: sequencing error, reverse-transcription misincorporation at modified RNA
bases, and NUMT-derived reads. That is the floor a real caller works against, so
that is what the model uses - position by position, because it is emphatically
not a constant (Phase 0 3.1). Two literature constants are carried alongside as
a sensitivity analysis, never instead of the measurement.

The error rate is estimated from the primary cohort only. The replication cohort
is not evidence for the map until Phase 8 (AGENTS.md 7.3).
"""
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
# One definition of the detection model, shared with the published tool. Two
# copies would drift, and this is exactly the code a subtle bug already hid in.
from mitofeasible import lod, min_alt_reads  # noqa: E402


def truth_nonref(vcf_path):
    """(sample, position) pairs where the donor's DNA carries a non-reference call.

    Everything absent from the callset is reference: that is what a sites VCF
    means, and it is also what makes the error estimate possible at all.
    """
    pairs, samples = set(), []
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("##"):
                continue
            p = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                samples = p[9:]
                continue
            pos = int(p[1])
            for s, g in zip(samples, p[9:]):
                gt = g.split(":")[0].replace("|", "/")
                if any(x not in ("0", ".") for x in gt.split("/")):
                    pairs.add((s, pos))
    return pairs, samples


def estimate_epsilon(alleles, samples, truth, cfg):
    """Per-position error rate, pooled across the primary cohort's donors.

    Only donor-positions whose DNA is reference contribute, so a real homoplasmic
    variant cannot be mistaken for noise. Pooling counts rather than averaging
    fractions keeps deep samples from being outvoted by shallow ones at the same
    position.
    """
    a = alleles[(alleles.arm == "A") & alleles.usable
                & (alleles.depth >= cfg["epsilon_min_depth"])]
    a = a.merge(samples[["run_accession", "sample_title", "bioproject"]],
                on="run_accession")
    a = a[a.bioproject == cfg["epsilon_source_cohort"]]
    keep = ~pd.MultiIndex.from_arrays([a.sample_title, a.pos]).isin(truth)
    a = a[keep]
    g = a.groupby("pos").agg(alt=("alt_count", "sum"), dep=("depth", "sum"),
                             n_donors=("run_accession", "nunique"),
                             max_donor_frac=("alt_frac", "max")).reset_index()
    g["epsilon"] = g.alt / g.dep
    return g


def build_map(alleles, samples, eps_tab, cfg):
    """Position x condition table of detectable allele fractions."""
    a = alleles[alleles.usable].merge(samples[["run_accession", "tissue"]],
                                      on="run_accession")
    cond = (a.groupby(["tissue", "arm", "pos"])
             .agg(median_depth=("depth", "median"), n_samples=("run_accession", "nunique"))
             .reset_index()
             .merge(eps_tab[["pos", "epsilon", "n_donors"]], on="pos", how="left"))
    alpha = cfg["alpha"] / cfg["n_tests"]
    # A position with no usable donor coverage gets the cohort's median rate
    # rather than a silently missing row - and is flagged so it can be excluded.
    fallback = float(eps_tab.epsilon.median())
    cond["epsilon_imputed"] = cond.epsilon.isna()
    cond["epsilon"] = cond.epsilon.fillna(fallback)
    cond["depth_used"] = np.floor(cond.median_depth).astype("int64")
    cond["min_detectable_af"] = lod(cond.depth_used.values, cond.epsilon.values,
                                    alpha, cfg["power"])
    for e in cfg["epsilon_sensitivity"]:
        cond[f"min_detectable_af_eps{e}"] = lod(cond.depth_used.values, e,
                                                alpha, cfg["power"])
    return cond


def selftest():
    # At 1% error and 1000x, a handful of alt reads is unremarkable; 40 is not.
    k = min_alt_reads(1000, 0.01, 1e-6)[0]
    assert 25 < k < 60, k
    # Deeper is never worse, and the limit must fall with depth.
    l100, l10000 = lod(100, 0.001, 1e-6, 0.8), lod(10000, 0.001, 1e-6, 0.8)
    assert l10000 < l100, (l100, l10000)
    assert 0 < l10000 < 0.05, l10000
    # A noisier position is harder, at the same depth.
    assert lod(1000, 0.01, 1e-6, 0.8) > lod(1000, 0.0001, 1e-6, 0.8)
    # Demanding more power can only raise the limit.
    assert lod(1000, 0.001, 1e-6, 0.95) >= lod(1000, 0.001, 1e-6, 0.8)
    assert np.isnan(lod(0, 0.001, 1e-6, 0.8)), "zero depth is not detectable"
    # A median depth over an even number of samples is a half-integer. scipy's
    # binom yields NaN for non-integer n without raising, which used to pin the
    # limit at 1.0 for every such position - half the map.
    half = lod(477.5, 4.68e-4, 1e-6, 0.8)
    assert np.isfinite(half) and half < 0.05, half
    assert abs(half - lod(477, 4.68e-4, 1e-6, 0.8)) < 1e-12, "half-integer must floor"

    import tempfile, os
    v = tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False)
    v.write("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA1\tNA2\n")
    v.write("MT\t73\t.\tA\tG\t.\t.\t.\tGT\t1\t0\n")
    v.write("MT\t263\t.\tA\tG\t.\t.\t.\tGT\t.\t1\n")
    v.close()
    pairs, names = truth_nonref(v.name)
    assert names == ["NA1", "NA2"]
    assert pairs == {("NA1", 73), ("NA2", 263)}, pairs
    os.unlink(v.name)

    al = pd.DataFrame({"run_accession": ["r1", "r1", "r2", "r2"], "arm": "A",
                       "usable": True, "pos": [73, 500, 73, 500],
                       "depth": [1000, 1000, 1000, 1000],
                       "alt_count": [990, 5, 3, 10],
                       "alt_frac": [.99, .005, .003, .01]})
    sm = pd.DataFrame({"run_accession": ["r1", "r2"], "sample_title": ["NA1", "NA2"],
                       "bioproject": ["P", "P"], "tissue": ["lcl", "lcl"]})
    cfg = {"epsilon_min_depth": 100, "epsilon_source_cohort": "P"}
    e = estimate_epsilon(al, sm, {("NA1", 73)}, cfg).set_index("pos")
    # r1 is a true variant at 73 and must be excluded; only r2's 3/1000 remains.
    assert abs(e.loc[73, "epsilon"] - 0.003) < 1e-9, e.loc[73].to_dict()
    assert abs(e.loc[500, "epsilon"] - 0.0075) < 1e-9, e.loc[500].to_dict()
    assert e.loc[73, "n_donors"] == 1
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/params.yaml")
    ap.add_argument("--samples", default="config/samples.tsv")
    ap.add_argument("--truth", default=None, help="subset VCF written by 06b")
    ap.add_argument("--out", default="results/tables")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    c = yaml.safe_load(open(a.config))
    f = c["feasibility"]
    samples = pd.read_csv(a.samples, sep="\t")
    files = sorted(Path(f["counts_dir"]).glob("*_alleles.parquet"))
    if not files:
        sys.exit(f"no allele counts in {f['counts_dir']}")
    print(f"reading {len(files)} allele tables", file=sys.stderr)
    alleles = pd.concat([pd.read_parquet(p) for p in files], ignore_index=True)

    truth_path = a.truth or Path(c["haplogroup"]["hsd_dir"]) / "truth_subset.vcf"
    truth, _ = truth_nonref(truth_path)
    eps = estimate_epsilon(alleles, samples, truth, f)
    cond = build_map(alleles, samples, eps, f)

    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    eps.to_csv(out / "error_rate_per_position.tsv", sep="\t", index=False)
    cond.to_parquet(out / "feasibility_map.parquet", index=False)
    q = eps.epsilon.quantile([.5, .9, .99, 1.0])
    print(f"epsilon: median {q[.5]:.2e}  p90 {q[.9]:.2e}  p99 {q[.99]:.2e}  max {q[1.0]:.2e}")
    print(cond.groupby(["tissue", "arm"]).min_detectable_af.median().to_string())


if __name__ == "__main__":
    main()
