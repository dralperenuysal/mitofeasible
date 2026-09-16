#!/usr/bin/env python3
"""Does library chemistry move the detection limit, and does it move it one way?

The chemistry cohort (PRJNA1086804, bulk DLPFC, poly(A) vs rRNA-depleted, the
same donors and centre either way) has no matched DNA, so its error floor is the
same RNA-only surrogate phase 8 validated in the primary cohort. Both arms of
the contrast are estimated the same way, so whatever bias the surrogate carries
is common to them and the comparison between them is the quantity of interest.

Nothing here reimplements the model: the estimator comes from phase 8 and the
detection limit from the released tool.
"""
import argparse, importlib.util, sys
from pathlib import Path
import numpy as np
import pandas as pd
import _config
from scipy import stats

HERE = Path(__file__).resolve().parent


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, HERE / path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def thin(alleles, target, seed=20260904):
    """Binomially thin each sample's counts to `target` reads per position.

    The surrogate is a median of per-donor non-reference fractions, so at low
    depth most donors observe zero alternate reads and the median collapses to
    exactly zero - the estimator reads as error-free where it is only
    undersampled. That bias depends on depth, and the two chemistries differ in
    depth roughly 18-fold, so the arms are not comparable until depth is matched.
    """
    rng = np.random.default_rng(seed)
    d = alleles.copy()
    keep = d.depth >= target
    d = d[keep].copy()
    frac = np.where(d.depth > 0, d.alt_count / d.depth, 0.0)
    d["alt_count"] = rng.binomial(target, np.clip(frac, 0, 1))
    d["depth"] = target
    d["alt_frac"] = d.alt_count / target
    return d


def pooled_epsilon(alleles, samples, cohort, cfg):
    """Pooled non-reference rate: total alternate reads over total depth.

    The median-of-donors surrogate needs each donor to resolve the floor on its
    own, which at a few hundred reads per donor it cannot: most donors observe no
    alternate read and the median is exactly zero. Pooling spends all donors'
    reads on one estimate, which is also how the DNA-anchored error rate in the
    primary cohort is computed. The median is still used, but only to drop fixed
    differences from the reference haplogroup.
    """
    a = alleles[(alleles.arm == "A") & alleles.usable
                & (alleles.depth >= cfg["surrogate_min_depth"])]
    a = a.merge(samples[["run_accession", "bioproject"]], on="run_accession")
    a = a[a.bioproject == cohort]
    g = a.groupby("pos").agg(alt=("alt_count", "sum"), dep=("depth", "sum"),
                             per_donor_median=("alt_frac", "median"),
                             n_donors=("run_accession", "nunique"),
                             median_depth=("depth", "median")).reset_index()
    g["surrogate"] = g.alt / g.dep
    g["fixed_difference"] = g.per_donor_median > cfg["surrogate_max_median_frac"]
    return g


def limits(alleles, samples, cfg, alpha, power, estimator="median"):
    """Per position, per chemistry: surrogate error floor and detectable AF."""
    rep = _load("rep08", "08_replication.py")
    mf = _load("mitofeasible", "mitofeasible.py")
    out = []
    for chem, g in samples.groupby("chemistry"):
        # surrogate_epsilon keys on a column named bioproject; the chemistry
        # label is what plays that role here, so it is passed under that name
        # rather than copying the estimator.
        s = g[["run_accession"]].assign(bioproject=chem)
        est = pooled_epsilon if estimator == "pooled" else rep.surrogate_epsilon
        e = est(alleles, s, chem, cfg)
        e = e[~e.fixed_difference].copy()
        e["min_detectable_af"] = mf.lod(e.median_depth.values, e.surrogate.values,
                                        alpha, power)
        out.append(e.assign(chemistry=chem))
    return pd.concat(out, ignore_index=True)


def yield_table(samples, stats_dir, arm="A"):
    """chrM share of alignments per sample, from the counts phase 3 wrote."""
    rows = []
    for r in samples.itertuples():
        f = Path(stats_dir) / f"{r.run_accession}_{arm}.counts"
        if not f.exists():
            continue
        d = pd.read_csv(f, sep="\t", names=["run", "arm", "what", "n"])
        n = dict(zip(d.what, d.n))
        rows.append({"run_accession": r.run_accession, "chemistry": r.chemistry,
                     "region": r.region, "total": n["total"], "chrM": n["chrM"],
                     "chrM_pct": n["chrM"] / n["total"] * 100})
    return pd.DataFrame(rows)


def paired_test(d, value, block="block", chem="chemistry"):
    """Wilcoxon signed-rank over blocks, because the design is paired.

    Each tissue block was sequenced both ways, so the two arms are not
    independent samples: donor, position and RNA all cancel within a block and
    only chemistry varies. An unpaired test here would throw that away and
    charge between-donor variance against the effect.
    """
    w = d.pivot_table(index=block, columns=chem, values=value)
    w = w.dropna()
    if len(w) < 3 or w.shape[1] != 2:
        return {"n_blocks": len(w)}
    a, b = w["polya"], w["rrna_depleted"]
    return {"n_blocks": len(w), "median_polya": float(a.median()),
            "median_rrna": float(b.median()),
            "median_paired_ratio": float((b / a).median()),
            "n_blocks_rrna_higher": int((b > a).sum()),
            "wilcoxon_p": float(stats.wilcoxon(a, b).pvalue)}


def by_gene_type(lim, genes):
    """Median limit per gene type, and the ratio between the two chemistries."""
    g = genes.copy()
    rows = []
    for r in g.itertuples():
        sel = lim[(lim.pos >= r.start) & (lim.pos <= r.end)]
        rows.append(sel.assign(gene=r.name, gene_type=r.gene_type))
    d = pd.concat(rows, ignore_index=True)
    t = d.pivot_table(index="gene_type", columns="chemistry",
                      values="min_detectable_af", aggfunc="median")
    if {"polya", "rrna_depleted"} <= set(t.columns):
        # Above 1: the rRNA-depleted library needs a higher allele fraction, so
        # it is the worse of the two at that class of positions.
        t["rrna_over_polya"] = t.rrna_depleted / t.polya
    return t.reset_index(), d


def selftest():
    mf = _load("mitofeasible", "mitofeasible.py")
    # Deeper coverage at the same error rate must not make detection harder.
    a = mf.lod(500.0, 1e-3, 1e-6, 0.8)
    b = mf.lod(5000.0, 1e-3, 1e-6, 0.8)
    assert b < a, (a, b)
    # A higher error floor at the same depth must.
    c = mf.lod(500.0, 1e-2, 1e-6, 0.8)
    assert c > a, (a, c)
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/params.yaml")
    ap.add_argument("--samples", default="config/samples_chem.tsv")
    ap.add_argument("--out", default="results/tables")
    # Pooled by default: this cohort's rRNA-depleted arm is too shallow per donor
    # for the median-of-donors estimator, which collapses to exactly zero at half
    # its positions. `median` is kept only to reproduce that diagnosis.
    ap.add_argument("--estimator", default="pooled", choices=["median", "pooled"],
                    help="how the DNA-free error floor is estimated")
    ap.add_argument("--downsample-to", type=int, default=None,
                    help="thin every sample to this per-position depth first")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()

    c = _config.load(a.config)
    f, r = c["feasibility"], c["replication"]
    alpha, power = f["alpha"] / f["n_tests"], f["power"]
    samples = pd.read_csv(a.samples, sep="\t")
    counts = Path(f["counts_dir"])
    frames = []
    for run in samples.run_accession:
        p = counts / f"{run}_alleles.parquet"
        if not p.exists():
            print(f"missing {p}", file=sys.stderr); continue
        frames.append(pd.read_parquet(p))
    alleles = pd.concat(frames, ignore_index=True)
    print(f"{alleles.run_accession.nunique()} samples", file=sys.stderr)

    tag = ""
    if a.downsample_to:
        n0 = len(alleles)
        alleles = thin(alleles, a.downsample_to)
        tag = f"_depth{a.downsample_to}"
        print(f"thinned to {a.downsample_to}x: {n0} -> {len(alleles)} rows",
              file=sys.stderr)

    tag += "" if a.estimator == "pooled" else "_median"
    lim = limits(alleles, samples, r, alpha, power, a.estimator)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    lim.to_csv(out / f"chemistry_limits_per_position{tag}.tsv", sep="\t", index=False)

    y = yield_table(samples, Path(c["alignment"]["bam_dir"]).parent / "stats")
    if len(y):
        # yield_table already carries region; only the donor is new here.
        y = y.merge(samples[["run_accession", "donor"]], on="run_accession")
        y["block"] = y.donor + "_" + y.region
        y.to_csv(out / "chemistry_yield.tsv", sep="\t", index=False)
        print("\nchrM %% of alignments:")
        print(y.groupby("chemistry").chrM_pct.agg(["count", "median", "min", "max"])
               .to_string(float_format="%.3f"))
        t = paired_test(y, "chrM_pct")
        print("paired over tissue blocks:", {k: (round(v, 6) if isinstance(v, float) else v)
                                             for k, v in t.items()})
        pd.DataFrame([t]).to_csv(out / "chemistry_yield_paired_test.tsv",
                                 sep="\t", index=False)

    genes = pd.read_csv(c["coverage"]["genes_bed"], sep="\t")
    t, per_gene = by_gene_type(lim, genes)
    t.to_csv(out / f"chemistry_limits_by_gene_type{tag}.tsv", sep="\t", index=False)
    print("\nmedian minimum detectable allele fraction:")
    print(t.to_string(index=False, float_format="%.5f"))
    print("\nsurrogate error floor (median) and depth:")
    print(lim.groupby("chemistry")[["surrogate", "median_depth"]].median()
             .to_string(float_format="%.6g"))
    print("\nzero-surrogate positions (the undersampling artefact):")
    print(lim.assign(z=lim.surrogate == 0).groupby("chemistry").z
             .agg(["size", "sum", "mean"]).to_string())
    print("\npositions supporting 1% detection:")
    print(lim.assign(ok=lim.min_detectable_af <= .01)
             .groupby("chemistry").ok.agg(["size", "sum", "mean"]).to_string())


if __name__ == "__main__":
    main()
