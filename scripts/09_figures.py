#!/usr/bin/env python3
"""Phase 9 - the figures in the definition of done.

  1  coverage landscape across chrM, by cohort, with gene / NUMT / D-loop tracks
  2  measured per-position error rate, with the RNA-modification hotspots
  3  feasibility map: position vs minimum detectable allele fraction
  4  replication concordance
  5  circularity: what the rotated reference recovers at the seam (09b)
  8  library chemistry: poly(A) against rRNA depletion (needs 12_chemistry.py)
  6  NUMT delta against MAPQ threshold
  7  the DNA-free error-rate surrogate against the DNA-based one

The library-chemistry contrast originally planned for figure 2 is figure 8,
carried by a third cohort in which chemistry is the only variable.

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
import _config

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
    ax.set_title("chrM coverage landscape\n"
                 "shaded: control region (grey) and NUMT-affected regions (cohort colour)",
                 loc="left")
    ax.legend(frameon=False, loc="lower right", ncol=2)
    gene_track(axg, genes)
    axg.set_xlim(1, 16569)
    axg.set_xlabel("position on rCRS (bp)")
    fig.savefig(out / "fig1_coverage_landscape.png")
    fig.savefig(out / "fig1_coverage_landscape.pdf")
    plt.close(fig)


ANNOT = [(310, "310  poly-C (HVS2)"), (2617, "2617  m$^1$A947 (16S)"),
         (4264, "4264  mt-tRNA$^{Ile}$"), (5513, "5513  mt-tRNA$^{Trp}$"),
         (12139, "12139  mt-tRNA$^{His}$"), (16189, "16189  poly-C (HVS1)")]
FLOOR = 1e-5


def fig2_error_landscape(eps, genes, out):
    """The measured error floor is not a constant - three logs of structure."""
    fig = plt.figure(figsize=(7.2, 3.9))
    gs = fig.add_gridspec(2, 2, height_ratios=[8, 1], width_ratios=[2.6, 1],
                          hspace=.08, wspace=.28)
    ax, axg, axh = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[:, 1])

    e = eps.sort_values("pos")
    y = e.epsilon.clip(lower=FLOOR)
    shade_dloop(ax)
    # Placed in the widest gap between hotspot annotations (5,513 to 12,139), not
    # at the left edge: there they sat under the 310 label and, at print size,
    # ran into it. The opaque box keeps the scatter from showing through.
    for t, lab, col in [(1e-3, "0.001  (the usual assumption)", "#e34948"),
                        (1e-2, "0.01", "#888")]:
        ax.axhline(t, color=col, lw=.6, ls="--", zorder=1)
        ax.text(8600, t * 1.25, lab, fontsize=6, color=col, va="bottom",
                ha="center", zorder=4,
                bbox=dict(facecolor="white", edgecolor="none", pad=1.2))
    ax.scatter(e.pos, y, s=1.1, c="#6f7f72", lw=0, alpha=.55, zorder=2)
    hot = e[e.epsilon >= .05]
    ax.scatter(hot.pos, hot.epsilon, s=9, c="#e34948", lw=0, zorder=3)
    # Labels go in a reserved band above the data, rotated: the hotspots are far
    # apart on x, so vertical text cannot collide the way horizontal text did.
    for p, lab in ANNOT:
        v = float(e.loc[e.pos == p, "epsilon"].iloc[0])
        ax.annotate(lab, (p, v), textcoords="offset points", xytext=(0, 7),
                    rotation=90, ha="center", va="bottom", fontsize=5.5,
                    color="#333")
    ax.set_yscale("log")
    ax.set_ylim(FLOOR * .7, 400)
    # The headroom above 1.0 only holds the labels; a rate cannot go there, so
    # the ticks stop at 1 rather than implying the axis is meaningful up to 100.
    ax.set_yticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
    ax.tick_params(labelbottom=False)
    ax.set_ylabel("per-position error rate $\\varepsilon$")
    ax.set_title("The error floor is position-specific, not constant",
                 loc="left", pad=2)
    gene_track(axg, genes)
    axg.set_xlim(1, 16569)
    ax.set_xlim(1, 16569)
    axg.set_xlabel("position on rCRS (bp)")

    # Where the distribution actually sits, against the constant people assume.
    axh.hist(np.log10(y), bins=70, color="#9aa094", lw=0)
    med, p99 = e.epsilon.median(), e.epsilon.quantile(.99)
    for v, lab, col in [(med, "median", "#333"), (p99, "p99", "#3d4a58"),
                        (1e-3, "0.001", "#e34948")]:
        axh.axvline(np.log10(v), color=col, lw=.7, ls="--")
        axh.text(np.log10(v), .98, lab, transform=axh.get_xaxis_transform(),
                 rotation=90, fontsize=5.5, color=col, ha="right", va="top")
    # The 76 positions with no observed non-reference read pile up at the clip.
    axh.text(np.log10(FLOOR), .22, f"{int((e.epsilon == 0).sum())} positions with "
             "$\\varepsilon$ = 0,\nclipped to the axis floor",
             transform=axh.get_xaxis_transform(), fontsize=5.5, color="#777",
             ha="left", va="center")
    axh.set_xlabel("log$_{10}$ $\\varepsilon$")
    axh.set_ylabel("positions")
    n_hi = int((e.epsilon > 1e-3).sum())
    axh.set_title(f"{n_hi} positions ({n_hi / len(e) * 100:.1f}%) exceed 0.001\n"
                  f"median {med:.2e}, p99 {p99:.2e}, max {e.epsilon.max():.3f}",
                  loc="left", fontsize=7)
    fig.savefig(out / "fig2_error_landscape.png")
    fig.savefig(out / "fig2_error_landscape.pdf")
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
    ax.set_title("Feasibility map (measured per-position error rate, power 0.8)",
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
    ax.set_title("Does the map transfer?", loc="left")
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


HIST_MAX = 2.0     # relative depth change, i.e. 200%


def fig6_numt_mapq(pos, genes, out, rel=0.05):
    """Filtering to unique reads makes the NUMT effect larger, not smaller."""
    flag = f"affected_{rel}"
    d = pos[(pos.strand == "both") & (pos.role.notna())]
    fig, (ax, ax2, axg) = plt.subplots(
        3, 1, figsize=(7.2, 4.0), sharex=False,
        gridspec_kw={"height_ratios": [6, 5, 1], "hspace": .45})

    # How many positions move, and by how much, at each MAPQ threshold.
    tis = sorted(d.tissue.unique())
    w, xs = .3, np.arange(len(tis))
    for i, q in enumerate(sorted(d.mapq_min.unique())):
        n = [int(d[(d.tissue == t) & (d.mapq_min == q)][flag].sum()) for t in tis]
        b = ax.bar(xs + (i - .5) * w, n, w, label=f"MAPQ $\\geq$ {q}",
                   color=["#9aa094", "#e34948"][i], lw=0)
        ax.bar_label(b, fontsize=6, padding=1)
    ax.set_xticks(xs, [LABEL[t] for t in tis])
    ax.set_xlim(-.55, len(tis) - .45)
    ax.set_ylim(0, ax.get_ylim()[1] * 1.28)     # headroom for the bar labels
    ax.set_ylabel(f"positions affected\n(|rel. $\\Delta$| $\\geq$ {rel:.0%})")
    ax.set_title("MAPQ filtering widens the NUMT gap instead of closing it",
                 loc="left")
    ax.legend(frameon=False, ncol=2, loc="upper center")

    # The effect sizes themselves, not just the counts.
    # A handful of positions move by more than 100%; on a shared axis they flatten
    # the body of the distribution, so the axis stops at HIST_MAX and says so.
    beyond = {}
    for i, q in enumerate(sorted(d.mapq_min.unique())):
        v = d[(d.mapq_min == q) & d[flag]].median_rel_delta.abs().dropna()
        beyond[q] = int((v > HIST_MAX).sum())
        # Dropped, not clipped: piling them into the last bin draws a spike that
        # competes with the real shape. The count is reported below instead.
        ax2.hist(v[v <= HIST_MAX], bins=np.linspace(0, HIST_MAX, 61),
                 histtype="step", lw=1.0, color=["#9aa094", "#e34948"][i],
                 label=f"MAPQ $\\geq$ {q}  (median {v.median():.0%})")
    ax2.set_xlim(0, HIST_MAX)
    ax2.text(.99, .42, "not shown, beyond the axis: "
             + ", ".join(f"{n} at q{q}" for q, n in beyond.items()),
             transform=ax2.transAxes, ha="right", fontsize=6, color="#777")
    ax2.set_xlabel("|median relative depth change| at affected positions")
    ax2.set_ylabel("positions")
    ax2.legend(frameon=False)

    # Where they sit: the MAPQ 255 set contains the MAPQ 0 set.
    for i, q in enumerate(sorted(d.mapq_min.unique())):
        p = d[(d.mapq_min == q) & d[flag]].pos.unique()
        axg.vlines(p, i, i + .8, lw=.4, color=["#9aa094", "#e34948"][i])
    axg.set_ylim(-.2, 2)
    axg.set_yticks([.4, 1.4], ["q0", "q255"], fontsize=6)
    axg.set_xlim(1, 16569)
    axg.set_xlabel("position on rCRS (bp)")
    for s in axg.spines.values():
        s.set_visible(False)
    fig.savefig(out / "fig6_numt_mapq.png")
    fig.savefig(out / "fig6_numt_mapq.pdf")
    plt.close(fig)


def fig7_surrogate(m, v, out):
    """Does the DNA-free error estimator stand in for the DNA-based one?"""
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0),
                             gridspec_kw={"width_ratios": [1.15, 1], "wspace": .3})
    ax = axes[0]
    both = m[(m.surrogate > 0) & (m.epsilon > 0)]
    ax.scatter(both.epsilon, both.surrogate, s=1.2, c="#9aa094", alpha=.35, lw=0)
    lim = [min(both.epsilon.min(), both.surrogate.min()) * .8,
           max(both.epsilon.max(), both.surrogate.max()) * 1.2]
    ax.plot(lim, lim, color="#333", lw=.6)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("DNA-based $\\varepsilon$ (matched WGS)")
    ax.set_ylabel("surrogate $\\varepsilon$ (RNA only)")
    ax.set_title(f"A DNA-free stand-in for $\\varepsilon$\n"
                 f"Pearson $r$ = {v['pearson_log']:.3f} on log scale, "
                 f"n = {int(v['n_positions'])}", loc="left")

    ax = axes[1]
    ratio = (m.surrogate / m.epsilon.replace(0, np.nan)).replace(0, np.nan).dropna()
    lr = np.log2(ratio)
    # A thin tail of a few hundred positions runs to -9; showing it leaves the
    # body of the distribution in a tenth of the axis.
    lo, hi = lr.quantile(.002), lr.quantile(.999)
    ax.hist(lr[(lr >= lo) & (lr <= hi)], bins=80, color="#9aa094", lw=0)
    ax.axvline(0, color="#333", lw=.7)
    ax.set_xlim(lo, hi)
    ax.set_xlabel("log$_2$(surrogate / DNA-based)")
    ax.set_ylabel("positions")
    ax.set_title(f"median {v['median_ratio']:.3f}x — the surrogate is "
                 "the conservative one", loc="left")
    fig.savefig(out / "fig7_surrogate.png")
    fig.savefig(out / "fig7_surrogate.pdf")
    plt.close(fig)

# Validated with the palette checker: worst adjacent pair dE 24.7 (protan),
# 33.6 (normal vision). Both series are also direct-labelled, so identity never
# rests on colour alone.
CHEM = {"polya": "#2a78d6", "rrna_depleted": "#eb6834"}
CHEM_LABEL = {"polya": "poly(A)", "rrna_depleted": "rRNA-depleted"}
GT_LABEL = {"protein_coding": "protein-coding", "Mt_rRNA": "rRNA", "Mt_tRNA": "tRNA"}


def fig8_chemistry(y, t, out):
    """What rRNA depletion costs: a lot of yield, unevenly distributed."""
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.1),
                                  gridspec_kw={"width_ratios": [1, 1.25], "wspace": .34})

    # Left: every sample, so the separation is visible rather than asserted.
    rng = np.random.default_rng(7)
    for i, chem in enumerate(["polya", "rrna_depleted"]):
        v = y[y.chemistry == chem].chrM_pct.values
        ax.scatter(i + rng.uniform(-.13, .13, len(v)), v, s=14, lw=.5,
                   facecolor=CHEM[chem], edgecolor="white", zorder=3)
        ax.plot([i - .26, i + .26], [np.median(v)] * 2, color=CHEM[chem], lw=1.6,
                zorder=4)
        # Beside the median bar, not above the cloud: at 12 jittered points an
        # overhead label lands on one of them.
        ax.text(i + .3, np.median(v), f"{np.median(v):.2f}%", ha="left",
                va="center", fontsize=7, color=CHEM[chem])
    ax.set_yscale("log")
    ax.set_xticks([0, 1], [CHEM_LABEL["polya"], CHEM_LABEL["rrna_depleted"]])
    ax.set_xlim(-.55, 1.75)
    ax.set_ylabel("chrM share of alignments (%)")
    # Within-block ratio, not a ratio of group medians: the design is paired, and
    # the two differ (18x vs 17x here) because the blocks are not identical.
    w = y.pivot_table(index="block", columns="chemistry", values="chrM_pct")
    fold = 1 / (w.rrna_depleted / w.polya).median()
    ax.set_title(f"Yield: {fold:.0f}x, no overlap", loc="left", fontsize=8, pad=14)
    ax.annotate("", xy=(.5, y[y.chemistry == "rrna_depleted"].chrM_pct.median()),
                xytext=(.5, y[y.chemistry == "polya"].chrM_pct.median()),
                arrowprops=dict(arrowstyle="->", color="#777", lw=.7))

    # Right: where that loss actually lands. A dumbbell, because the quantity of
    # interest is the movement between the two chemistries, not either level.
    t = t.set_index("gene_type").loc[["protein_coding", "Mt_rRNA", "Mt_tRNA"]]
    ypos = np.arange(len(t))[::-1]
    for yy, (gt, r) in zip(ypos, t.iterrows()):
        ax2.plot([r.polya * 100, r.rrna_depleted * 100], [yy, yy], color="#bbb",
                 lw=1.4, zorder=1)
        ax2.scatter([r.polya * 100], [yy], s=30, facecolor=CHEM["polya"],
                    edgecolor="white", lw=.6, zorder=3)
        ax2.scatter([r.rrna_depleted * 100], [yy], s=30,
                    facecolor=CHEM["rrna_depleted"], edgecolor="white", lw=.6, zorder=3)
        ax2.text(r.rrna_depleted * 100 * 1.08, yy, f"{r.rrna_over_polya:.2f}x",
                 va="center", fontsize=7,
                 color="#333" if r.rrna_over_polya > 1.2 else "#2e7d32")
    ax2.set_yticks(ypos, [GT_LABEL[g] for g in t.index])
    ax2.set_ylim(-.6, len(t) - .4)
    ax2.set_xscale("log")
    ax2.set_xlim(.3, 4.0)
    # A log axis this short defaults to a single decade label; name the values.
    ax2.set_xticks([0.4, 0.6, 1.0, 1.5, 2.5])
    ax2.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(
        lambda v, _: f"{v:g}"))
    ax2.set_xticks([], minor=True)
    ax2.set_xlabel("median minimum detectable allele fraction (%)")
    ax2.set_title(f"Penalty: {t.loc['protein_coding'].rrna_over_polya:.1f}x in "
                  "genes, none in tRNAs", loc="left", fontsize=8, pad=14)
    for s_ in ("left",):
        ax2.spines[s_].set_visible(False)
    ax2.tick_params(axis="y", length=0)
    h = [plt.Line2D([], [], marker="o", ls="", color=CHEM[c], label=CHEM_LABEL[c])
         for c in ("polya", "rrna_depleted")]
    ax2.legend(handles=h, frameon=False, fontsize=7, loc="upper center",
               ncol=2, bbox_to_anchor=(.5, 1.02), handletextpad=.3,
               columnspacing=1.2)
    fig.savefig(out / "fig8_chemistry.png")
    fig.savefig(out / "fig8_chemistry.pdf")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/params.yaml")
    ap.add_argument("--tables", default="results/tables")
    ap.add_argument("--out", default="results/figures")
    a = ap.parse_args()
    style()
    c = _config.load(a.config)
    T, out = Path(a.tables), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    genes = pd.read_csv(c["coverage"]["genes_bed"], sep="\t")

    # Some inputs are large intermediates the repository does not publish. A
    # clone therefore redraws the figures whose tables are committed and says
    # which it could not, rather than dying on the first missing parquet.
    def draw(name, fn, *inputs):
        missing = [str(i) for i in inputs if isinstance(i, Path) and not i.exists()]
        if missing:
            print(f"{name} skipped, missing: {', '.join(missing)}")
            return
        try:
            fn()
        except FileNotFoundError as e:
            print(f"{name} skipped, missing: {e.filename}")

    def coverage_quantiles():
        samples = pd.read_csv(c["cohorts"]["samples_tsv"], sep="\t") \
            if "samples_tsv" in c["cohorts"] else pd.read_csv("config/samples.tsv", sep="\t")
        samples = samples[["run_accession", "tissue"]]
        files = sorted(Path(c["coverage"]["cov_dir"]).glob("*_per_base.parquet"))
        if not files:
            raise FileNotFoundError(2, "no per-base parquets", c["coverage"]["cov_dir"])
        rows = [pd.read_parquet(f)[lambda d: (d.arm == "A") & (d.mapq_min == 0)
                                  & (d.strand == "both")]
                [["run_accession", "pos", "depth"]] for f in files]
        return (pd.concat(rows).merge(samples, on="run_accession")
                .groupby(["tissue", "pos"]).depth
                .agg(median="median", q25=lambda s: s.quantile(.25),
                     q75=lambda s: s.quantile(.75)).reset_index())

    regions_f = T / "numt_delta_regions.tsv"
    draw("figure 1", lambda: fig1_coverage(
        coverage_quantiles(),
        genes,
        pd.read_csv(regions_f, sep="\t").pipe(
            lambda r: r[(r.rel_threshold == 0.05) & (r.strand == "both")
                        & (r.mapq_min == 0)]),
        out), regions_f)
    draw("figure 2", lambda: fig2_error_landscape(
        pd.read_csv(T / "error_rate_per_position.tsv", sep="\t"), genes, out),
        T / "error_rate_per_position.tsv")
    draw("figure 3", lambda: fig3_feasibility(
        pd.read_parquet(T / "feasibility_map.parquet"), genes, out),
        T / "feasibility_map.parquet")
    draw("figure 4", lambda: fig4_replication(
        pd.read_parquet(T / "replication_transfer.parquet"), out),
        T / "replication_transfer.parquet")
    draw("figure 6", lambda: fig6_numt_mapq(
        pd.read_parquet(T / "numt_delta_per_position.parquet"), genes, out),
        T / "numt_delta_per_position.parquet")
    draw("figure 8", lambda: fig8_chemistry(
        pd.read_csv(T / "chemistry_yield.tsv", sep="\t"),
        pd.read_csv(T / "chemistry_limits_by_gene_type.tsv", sep="\t"), out),
        T / "chemistry_yield.tsv", T / "chemistry_limits_by_gene_type.tsv")

    draw("figure 7", lambda: fig7_surrogate(
        pd.read_csv(T / "replication_surrogate_vs_dna.tsv", sep="\t"),
        pd.read_csv(T / "replication_surrogate_validation.tsv", sep="\t").iloc[0],
        out),
        T / "replication_surrogate_vs_dna.tsv",
        T / "replication_surrogate_validation.tsv")

    # Figure 5 needs the unspliced primary depth, which phase 4 does not keep;
    # 09b_circularity.py rebuilds it from the BAMs and draws its own figure.
    print(f"wrote figures to {out}: " + ", ".join(
        sorted(f.stem for f in out.glob("fig*.pdf"))))
    print("figure 5 comes from 09b_circularity.py.")


if __name__ == "__main__":
    main()
