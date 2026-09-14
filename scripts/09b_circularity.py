#!/usr/bin/env python3
"""Figure 5 - what the rotated reference recovers at the 16569/1 seam.

Phase 4 splices the control region in from the shifted pass and keeps only the
result, so the unspliced primary depth needed for a before/after picture is
gone by the time the tables are written. This rebuilds both profiles straight
from the chrM BAMs (they are small; only chrM is read), averages them per
cohort, and draws the figure. Nothing here feeds the feasibility map - it
documents an alignment decision made in phase 3.
"""
import argparse, importlib.util, sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent
C = {"lcl": "#2a78d6", "skeletal_muscle": "#eb6834"}
LABEL = {"lcl": "LCL (n=30)", "skeletal_muscle": "Skeletal muscle (n=30)"}
SEAM_FLANK = 250           # bp either side of the seam to show
SEAM_WINDOW = 10           # bp either side used for the headline depth ratio


def _cov():
    spec = importlib.util.spec_from_file_location("cov04", HERE / "04_coverage.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def seam_coords(length, flank):
    """Positions around the seam, as (1-based position, signed seam offset).

    Offset 0 is position `length`; +1 is position 1. Plotting against the offset
    puts the junction in the middle of the axis instead of at both edges.
    """
    left = np.arange(length - flank + 1, length + 1)      # ... 16569
    right = np.arange(1, flank + 1)                       # 1 ...
    pos = np.concatenate([left, right])
    off = np.concatenate([left - length, right])
    return pos, off


def profile(cfg, ref, aln, samples, arm="A", q=0):
    """Mean primary and shifted depth per position, per cohort."""
    cov = _cov()
    length, shift = ref["length"], ref["shift"]
    bam_dir = Path(aln["bam_dir"])
    rows = []
    for r in samples.itertuples():
        prim = cov.bam_depth(bam_dir / f"{r.run_accession}_{arm}_chrM.bam",
                             length, [q], False)[(q, "both")]
        shft = cov.bam_depth(
            bam_dir / f"{r.run_accession}_{arm}_shifted_Aligned.sortedByCoord.out.bam",
            length, [q], False, offset=shift)[(q, "both")]
        rows.append(pd.DataFrame({"tissue": r.tissue, "pos": np.arange(1, length + 1),
                                  "primary": prim, "shifted": shft}))
        print(f"  {r.run_accession}", file=sys.stderr)
    return (pd.concat(rows).groupby(["tissue", "pos"])[["primary", "shifted"]]
            .mean().reset_index())


def fig5_circularity(p, length, window, out):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0),
                             gridspec_kw={"width_ratios": [1.3, 1], "wspace": .3})
    pos, off = seam_coords(length, SEAM_FLANK)
    order = pd.DataFrame({"pos": pos, "off": off})

    ax = axes[0]
    ax.axvline(0.5, color="#e34948", lw=.7, ls="--")
    ax.text(0.5, .985, "16569 / 1", color="#e34948", fontsize=6, ha="center",
            va="top", transform=ax.get_xaxis_transform())
    for tis, d in p.groupby("tissue"):
        d = order.merge(d, on="pos").sort_values("off")
        ax.plot(d.off, d.primary, color=C[tis], lw=.9, ls=":")
        ax.plot(d.off, d.shifted, color=C[tis], lw=.9,
                label=LABEL[tis].split(" (")[0])
        ax.fill_between(d.off, d.primary, d.shifted, where=d.shifted > d.primary,
                        color=C[tis], alpha=.13, lw=0)
    ax.set_yscale("log")
    ax.set_xlabel("distance from the seam (bp)")
    ax.set_ylabel("mean depth")
    # Two encodings, two keys: colour is the cohort, line style is the pass.
    style_keys = [plt.Line2D([], [], color="#555", lw=.9, ls=":", label="linear rCRS"),
                  plt.Line2D([], [], color="#555", lw=.9, label="rotated 8 kb")]
    ax.add_artist(ax.legend(handles=style_keys, frameon=False, fontsize=6,
                            loc="lower left"))
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    g = seam_gain(p, length)
    ax.set_title("The rotated pass recovers the seam\n"
                 + ", ".join(f"{LABEL[t].split(' (')[0]} {r.gain:.2f}x over "
                             f"$\\pm${SEAM_WINDOW} bp"
                             for t, r in g.set_index("tissue").iterrows()),
                 loc="left", fontsize=8)

    ax = axes[1]
    lo, hi = window
    d = p[(p.pos >= lo) & (p.pos <= hi)]
    for tis, g in d.groupby("tissue"):
        ax.scatter(g.primary, g.shifted, s=1.2, c=C[tis], lw=0, alpha=.4)
    lim = [max(d[["primary", "shifted"]].min().min(), .5),
           d[["primary", "shifted"]].max().max() * 1.2]
    ax.plot(lim, lim, color="#333", lw=.6)
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lim); ax.set_ylim(lim)
    r = d.groupby("tissue").apply(
        lambda g: g.primary.corr(g.shifted, method="spearman"))
    ax.set_xlabel("linear rCRS depth")
    ax.set_ylabel("rotated reference depth")
    ax.set_title(f"Seam-free window {lo}-{hi} bp\n"
                 + ", ".join(f"{LABEL[t].split(' (')[0]} $\\rho$ = {v:.5f}"
                             for t, v in r.items()), loc="left", fontsize=7)
    fig.savefig(out / "fig5_circularity.png")
    fig.savefig(out / "fig5_circularity.pdf")
    plt.close(fig)


def seam_gain(p, length, w=SEAM_WINDOW):
    """Mean depth over the w bp either side of the seam, both passes."""
    sel = (p.pos > length - w) | (p.pos <= w)
    g = p[sel].groupby("tissue")[["primary", "shifted"]].mean()
    g["gain"] = g.shifted / g.primary
    return g.reset_index()


def selftest():
    pos, off = seam_coords(16569, 5)
    assert list(pos) == [16565, 16566, 16567, 16568, 16569, 1, 2, 3, 4, 5]
    assert list(off) == [-4, -3, -2, -1, 0, 1, 2, 3, 4, 5]
    assert (np.diff(off) == 1).all(), "the seam axis must be contiguous"
    p = pd.DataFrame({"tissue": "lcl", "pos": np.arange(1, 16570),
                      "primary": 100.0, "shifted": 300.0})
    assert abs(seam_gain(p, 16569).gain.iloc[0] - 3.0) < 1e-9
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/params.yaml")
    ap.add_argument("--samples", default="config/samples.tsv")
    ap.add_argument("--tables", default="results/tables")
    ap.add_argument("--out", default="results/figures")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()

    c = yaml.safe_load(open(a.config))
    ref, aln, cov = c["reference"], c["alignment"], c["coverage"]
    T, out = Path(a.tables), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    cache = T / "shift_profile.tsv"
    if cache.exists():
        p = pd.read_csv(cache, sep="\t")
        print(f"reusing {cache}", file=sys.stderr)
    else:
        samples = pd.read_csv(a.samples, sep="\t")
        p = profile(c["coverage"], ref, aln, samples)
        p.to_csv(cache, sep="\t", index=False)

    plt.rcParams.update({
        "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 9,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.linewidth": .6, "figure.dpi": 150, "savefig.dpi": 300,
        "savefig.bbox": "tight"})
    fig5_circularity(p, ref["length"], cov["concordance_window"], out)
    g = seam_gain(p, ref["length"])
    g.to_csv(T / "shift_seam_gain.tsv", sep="\t", index=False)
    print(g.to_string(index=False))
    print(f"wrote figure 5 to {out}")


if __name__ == "__main__":
    main()
