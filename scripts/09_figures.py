#!/usr/bin/env python3
"""Phase 9 - the figures in the definition of done.

  1  coverage landscape across chrM, by cohort, with gene / NUMT / D-loop tracks
  2  library chemistry contrast - not produced: both cohorts are poly(A)
  3  feasibility map: position vs minimum detectable allele fraction
  4  replication concordance

Everything is drawn from the tables the earlier phases wrote; nothing is
recomputed here, so a figure can never disagree with the table behind it.
"""
import argparse, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

C = {"lcl": "#2a78d6", "skeletal_muscle": "#eb6834"}
LABEL = {"lcl": "LCL (n=30)", "skeletal_muscle": "Skeletal muscle (n=30)"}
GT = {"Mt_rRNA": "#3d4a58", "protein_coding": "#6f7f72", "Mt_tRNA": "#b9c2b4"}
DLOOP = [(16024, 16569), (1, 576)]


def style():
    plt.rcParams.update({
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": .6, "xtick.major.width": .6, "ytick.major.width": .6,
        "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
    })


def gene_track(ax, genes, y=0, h=1):
    """Gene rectangles, strand-split, on their own thin axis."""
    for g in genes.itertuples():
        ax.add_patch(plt.Rectangle((g.start, y if g.strand == "+" else y + h / 2),
                                   g.end - g.start, h / 2 * .85,
                                   facecolor=GT[g.gene_type], edgecolor="none"))
        if g.end - g.start > 700:
            ax.text((g.start + g.end) / 2, y - .25, g.name.replace("MT-", ""),
                    ha="center", va="top", fontsize=5.5, color="#444")
    ax.set_ylim(-1.1, h)
    ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


def shade_dloop(ax):
    for s, e in DLOOP:
        ax.axvspan(s, e, color="#000", alpha=.045, lw=0, zorder=0)


def fig1_coverage(cov, genes, regions, out):
    """Median depth per position with the interquartile band, per cohort."""
    fig, (ax, axg) = plt.subplots(2, 1, figsize=(7.2, 3.4), sharex=True,
                                  gridspec_kw={"height_ratios": [8, 1], "hspace": .08})
    shade_dloop(ax)
    for r in regions.itertuples():
        ax.axvspan(r.start, r.end, color=C[r.tissue], alpha=.07, lw=0, zorder=0)
    for tis, d in cov.groupby("tissue"):
        d = d.sort_values("pos")
        # d["median"], not d.median: the latter is pandas' method, and the
        # attribute form silently hands matplotlib a bound method.
        ax.fill_between(d.pos, d["q25"], d["q75"], color=C[tis], alpha=.22, lw=0)
        ax.plot(d.pos, d["median"], color=C[tis], lw=.7, label=LABEL[tis])
    ax.set_yscale("log")
    ax.set_ylabel("depth (x)")
    ax.set_title("Figure 1  chrM coverage landscape\n"
                 "shaded: control region (grey) and NUMT-affected regions (cohort colour)",
                 loc="left")
    ax.legend(frameon=False, loc="lower right", ncol=2)
    gene_track(axg, genes)
    axg.set_xlim(1, 16569)
    axg.set_xlabel("position on rCRS (bp)")
    fig.savefig(out / "fig1_coverage_landscape.png")
    fig.savefig(out / "fig1_coverage_landscape.pdf")
    plt.close(fig)


def fig3_feasibility(fmap, genes, out):
    """The core deliverable: what each position can detect."""
    fig, (ax, axg) = plt.subplots(2, 1, figsize=(7.2, 3.4), sharex=True,
                                  gridspec_kw={"height_ratios": [8, 1], "hspace": .08})
    shade_dloop(ax)
    for t, lab in [(.01, "1%"), (.05, "5%")]:
        ax.axhline(t, color="#888", lw=.5, ls=":", zorder=1)
        ax.text(16500, t, lab, fontsize=6, color="#777", va="bottom", ha="right")
    for tis, d in fmap[fmap.arm == "A"].groupby("tissue"):
        d = d.sort_values("pos")
        ax.plot(d.pos, d.min_detectable_af, color=C[tis], lw=.6, label=LABEL[tis])
    ax.set_yscale("log")
    ax.set_ylabel("min detectable allele fraction")
    ax.set_title("Figure 3  Feasibility map (measured per-position error rate, power 0.8)",
                 loc="left")
    ax.legend(frameon=False, loc="upper right", ncol=2)
    gene_track(axg, genes)
    axg.set_xlim(1, 16569)
    axg.set_xlabel("position on rCRS (bp)")
    fig.savefig(out / "fig3_feasibility_map.png")
    fig.savefig(out / "fig3_feasibility_map.pdf")
    plt.close(fig)


def fig4_replication(t, out):
    """Predicted vs observed limits on the held-out cohort, and the failures."""
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0),
                             gridspec_kw={"width_ratios": [1.15, 1], "wspace": .3})
    ax = axes[0]
    ok = t[~t.optimistic]
    ax.scatter(ok.af_predicted, ok.af_observed, s=1.2, c="#9aa094", alpha=.35, lw=0,
               label="transfers")
    bad = t[t.optimistic]
    ax.scatter(bad.af_predicted, bad.af_observed, s=4, c="#e34948", lw=0,
               label=f"optimistic ({len(bad)})")
    lim = [min(t.af_predicted.min(), t.af_observed.min()) * .8,
           max(t.af_predicted.max(), t.af_observed.max()) * 1.2]
    ax.plot(lim, lim, color="#333", lw=.6)
    ax.plot(lim, [x * 2 for x in lim], color="#333", lw=.5, ls="--")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("predicted from the primary cohort")
    ax.set_ylabel("observed in the replication cohort")
    ax.set_title("Figure 4  Does the map transfer?", loc="left")
    ax.legend(frameon=False, loc="upper left")

    ax = axes[1]
    r = np.log2(t.ratio.replace(0, np.nan).dropna())
    ax.hist(r, bins=90, color="#9aa094", lw=0)
    ax.axvline(0, color="#333", lw=.7)
    ax.axvline(1, color="#e34948", lw=.7, ls="--")
    ax.set_xlabel("log2(observed / predicted)")
    ax.set_ylabel("positions")
    ax.set_title(f"median {t.ratio.median():.2f}x — conservative side", loc="left")
    fig.savefig(out / "fig4_replication.png")
    fig.savefig(out / "fig4_replication.pdf")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/params.yaml")
    ap.add_argument("--tables", default="results/tables")
    ap.add_argument("--out", default="results/figures")
    a = ap.parse_args()
    style()
    c = yaml.safe_load(open(a.config))
    T, out = Path(a.tables), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(c["coverage"]["genes_bed"], sep="\t")

    # Coverage quantiles across each cohort, from the per-sample parquets.
    samples = pd.read_csv("config/samples.tsv", sep="\t")[["run_accession", "tissue"]]
    rows = []
    for f in sorted(Path(c["coverage"]["cov_dir"]).glob("*_per_base.parquet")):
        d = pd.read_parquet(f)
        rows.append(d[(d.arm == "A") & (d.mapq_min == 0) & (d.strand == "both")]
                    [["run_accession", "pos", "depth"]])
    cov = (pd.concat(rows).merge(samples, on="run_accession")
           .groupby(["tissue", "pos"]).depth
           .agg(median="median", q25=lambda s: s.quantile(.25),
                q75=lambda s: s.quantile(.75)).reset_index())

    regions = pd.read_csv(T / "numt_delta_regions.tsv", sep="\t")
    regions = regions[(regions.rel_threshold == 0.05) & (regions.strand == "both")
                      & (regions.mapq_min == 0)]
    fig1_coverage(cov, genes, regions, out)
    fig3_feasibility(pd.read_parquet(T / "feasibility_map.parquet"), genes, out)
    fig4_replication(pd.read_parquet(T / "replication_transfer.parquet"), out)
    print(f"wrote figures 1, 3 and 4 to {out}")
    print("figure 2 (library chemistry contrast) is not produced: both cohorts "
          "are poly(A); see README.")


if __name__ == "__main__":
    main()
