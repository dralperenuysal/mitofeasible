#!/usr/bin/env python3
"""Phase 8 - does the feasibility map transfer to an independent cohort?

The map has two halves. Depth is read from whatever data a user has and does not
transfer by construction. The per-position error rate is the part this study
claims is transferable, so that is what is tested here.

The muscle donors have no matched DNA, so their floor cannot be measured the way
GEUVADIS's was. It is estimated instead from the median non-reference fraction
across 30 unrelated donors: a genuine variant is carried by a minority, so the
median is dominated by error. That surrogate is first checked against the
DNA-based estimate in the primary cohort - if it cannot reproduce a number we
already know, it cannot be trusted on a cohort where we do not.

A partial failure is a result. Nothing here is tuned to improve agreement.
"""
import argparse, sys
from pathlib import Path
import numpy as np
import pandas as pd
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from mitofeasible import lod  # noqa: E402


def surrogate_epsilon(alleles, samples, cohort, cfg):
    """DNA-free error floor: the median non-reference fraction across donors.

    Positions whose median sits above `surrogate_max_median_frac` are fixed
    differences from rCRS - the reference is one particular haplogroup, so most
    people carry a few dozen - and are excluded rather than counted as noise.
    """
    a = alleles[(alleles.arm == "A") & alleles.usable
                & (alleles.depth >= cfg["surrogate_min_depth"])]
    a = a.merge(samples[["run_accession", "bioproject"]], on="run_accession")
    a = a[a.bioproject == cohort]
    g = a.groupby("pos").agg(surrogate=("alt_frac", "median"),
                             n_donors=("run_accession", "nunique"),
                             median_depth=("depth", "median")).reset_index()
    g["fixed_difference"] = g.surrogate > cfg["surrogate_max_median_frac"]
    return g


def validate_surrogate(surr_primary, eps_dna):
    """Does the DNA-free estimator reproduce the DNA-based one, where both exist?"""
    m = surr_primary.merge(eps_dna[["pos", "epsilon"]], on="pos", how="inner")
    m = m[~m.fixed_difference]
    both = m[(m.surrogate > 0) & (m.epsilon > 0)]
    return {
        "n_positions": int(len(m)),
        "spearman": float(m.surrogate.corr(m.epsilon, method="spearman")),
        "pearson_log": float(np.log10(both.surrogate).corr(np.log10(both.epsilon))),
        "median_surrogate": float(m.surrogate.median()),
        "median_dna_epsilon": float(m.epsilon.median()),
        "median_ratio": float((m.surrogate / m.epsilon.replace(0, np.nan)).median()),
    }, m


def transfer(eps_primary, surr_repl, cfg, alpha, power):
    """Predicted vs observed limits on the replication cohort's own depth.

    Depth is held identical between the two - it is the replication cohort's
    measured depth in both - so any difference is the error rate transferring or
    failing to, and nothing else.
    """
    m = surr_repl[~surr_repl.fixed_difference].merge(
        eps_primary[["pos", "epsilon"]], on="pos", how="inner")
    d = m.median_depth.values
    m["af_predicted"] = lod(d, m.epsilon.values, alpha, power)
    # The surrogate is a fraction, not a rate estimated from pooled counts; used
    # as eps it is the cohort's own floor at that position.
    m["af_observed"] = lod(d, np.maximum(m.surrogate.values, 1e-6), alpha, power)
    m["ratio"] = m.af_observed / m.af_predicted
    m["optimistic"] = m.ratio >= cfg["optimism_factor"]
    m["conservative"] = m.ratio <= 1 / cfg["optimism_factor"]
    return m


def selftest():
    al = pd.DataFrame({
        "run_accession": ["r1", "r2", "r3"] * 3, "arm": "A", "usable": True,
        "pos": [10] * 3 + [20] * 3 + [30] * 3,
        "depth": [1000] * 9,
        "alt_frac": [0.001, 0.002, 0.003,      # pos 10: noise
                     0.99, 0.98, 0.99,          # pos 20: fixed difference vs rCRS
                     0.001, 0.90, 0.002]})      # pos 30: one donor carries a variant
    sm = pd.DataFrame({"run_accession": ["r1", "r2", "r3"], "bioproject": ["P"] * 3})
    cfg = {"surrogate_min_depth": 100, "surrogate_max_median_frac": 0.10,
           "optimism_factor": 2.0}
    s = surrogate_epsilon(al, sm, "P", cfg).set_index("pos")
    assert abs(s.loc[10, "surrogate"] - 0.002) < 1e-9
    assert s.loc[20, "fixed_difference"], "a near-universal difference is not noise"
    # One carrier out of three must not drag the floor up: that is the point of
    # using the median rather than the mean.
    assert abs(s.loc[30, "surrogate"] - 0.002) < 1e-9, s.loc[30].to_dict()
    assert not s.loc[30, "fixed_difference"]

    eps = pd.DataFrame({"pos": [10, 30], "epsilon": [0.002, 0.0001]})
    t = transfer(eps, s.reset_index(), cfg, 1e-6, 0.8).set_index("pos")
    # Same epsilon on both sides -> the limits must coincide exactly.
    assert abs(t.loc[10, "ratio"] - 1.0) < 1e-6, t.loc[10].to_dict()
    # Primary says the position is 20x cleaner than it is -> optimistic.
    assert t.loc[30, "optimistic"] and not t.loc[30, "conservative"], t.loc[30].to_dict()

    v, _ = validate_surrogate(s.reset_index(), pd.DataFrame({"pos": [10, 30],
                                                             "epsilon": [0.002, 0.002]}))
    assert v["n_positions"] == 2 and abs(v["median_ratio"] - 1.0) < 1e-9, v
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
    c = yaml.safe_load(open(a.config))
    f, r = c["feasibility"], c["replication"]
    samples = pd.read_csv(a.samples, sep="\t")
    files = sorted(Path(f["counts_dir"]).glob("*_alleles.parquet"))
    print(f"reading {len(files)} allele tables", file=sys.stderr)
    alleles = pd.concat([pd.read_parquet(p) for p in files], ignore_index=True)

    primary = f["epsilon_source_cohort"]
    repl = next(k for k in c["cohorts"]["studies"] if k != primary)
    eps_dna = pd.read_csv(Path(a.out) / "error_rate_per_position.tsv", sep="\t")

    surr_p = surrogate_epsilon(alleles, samples, primary, r)
    v, merged = validate_surrogate(surr_p, eps_dna)
    print("surrogate vs DNA-based epsilon in the primary cohort:")
    for k, x in v.items():
        print(f"  {k}: {x:.4g}" if isinstance(x, float) else f"  {k}: {x}")

    surr_r = surrogate_epsilon(alleles, samples, repl, r)
    alpha = f["alpha"] / f["n_tests"]
    t = transfer(eps_dna, surr_r, r, alpha, f["power"])

    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([v]).to_csv(out / "replication_surrogate_validation.tsv",
                             sep="\t", index=False)
    merged.to_csv(out / "replication_surrogate_vs_dna.tsv", sep="\t", index=False)
    t.to_parquet(out / "replication_transfer.parquet", index=False)
    # The optimistic positions are the failure mode the paper asks readers to
    # act on, so they ship as a plain table rather than only inside a parquet
    # that the repository does not publish.
    cols = ["pos", "median_depth", "epsilon", "surrogate", "af_predicted",
            "af_observed", "ratio"]
    (t[t.optimistic].sort_values("ratio", ascending=False)[cols]
     .to_csv(out / "replication_optimistic_positions.tsv", sep="\t", index=False))
    print(f"\ntransfer on {len(t)} positions of the replication cohort:")
    print(f"  spearman(predicted, observed) = {t.af_predicted.corr(t.af_observed, method='spearman'):.4f}")
    print(f"  median ratio observed/predicted = {t.ratio.median():.3f}")
    print(f"  optimistic (>= {r['optimism_factor']}x too low): {int(t.optimistic.sum())} "
          f"({t.optimistic.mean()*100:.2f}%)")
    print(f"  conservative: {int(t.conservative.sum())} ({t.conservative.mean()*100:.2f}%)")


if __name__ == "__main__":
    main()
