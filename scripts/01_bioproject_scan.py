#!/usr/bin/env python3
"""Phase 1 - scan ENA for bulk RNA-seq BioProjects meeting the design constraints.

Downloads nothing but metadata. Writes results/tables/bioproject_candidates.tsv
and stops. Cohort selection is a human decision (AGENTS.md 7.2).
"""
import argparse, sys, urllib.parse, urllib.request
from pathlib import Path

import pandas as pd
import _config

FIELDS = [
    "study_accession", "sample_accession", "run_accession", "instrument_model",
    "library_layout", "library_selection", "read_count", "base_count",
    "study_title", "sample_title", "sample_description", "tissue_type",
    "description", "first_public",
]

# Tissue is recorded inconsistently at ENA: sample_description carries it ~20x
# more often than tissue_type (691 vs 34 runs for "heart"), so search all three.
TISSUE_FIELDS = ["sample_title", "sample_description", "tissue_type"]


def ena_query(cfg, selection, tissue):
    """One ENA read_run query for a (library_selection, tissue) cell."""
    q = (
        f'tax_id={cfg["tax_id"]}'
        f' AND library_strategy="{cfg["library_strategy"]}"'
        f' AND library_source="{cfg["library_source"]}"'
        f' AND instrument_platform="{cfg["instrument_platform"]}"'
        f' AND library_selection="{selection}"'
        " AND (" + " OR ".join(f'{f}="*{tissue}*"' for f in TISSUE_FIELDS) + ")"
    )
    url = cfg["ena_api"] + "?" + urllib.parse.urlencode({
        "result": "read_run", "query": q, "fields": ",".join(FIELDS),
        "format": "tsv", "limit": cfg["max_runs_fetched"],
    })
    with urllib.request.urlopen(url, timeout=300) as r:
        text = r.read().decode()
    if not text.strip():
        return pd.DataFrame(columns=FIELDS)
    df = pd.read_csv(pd.io.common.StringIO(text), sep="\t", dtype=str)
    df["library_selection"] = selection
    df["tissue"] = tissue
    return df


def annotate(df, cfg):
    """Derive the columns the design decisions are actually made on."""
    df = df.copy()
    df["read_count"] = pd.to_numeric(df["read_count"], errors="coerce")
    df["base_count"] = pd.to_numeric(df["base_count"], errors="coerce")
    # ENA gives no read-length field; derive it. Paired runs report both mates.
    df["read_length"] = (df["base_count"] / df["read_count"].where(df["read_count"] > 0)).round()
    df.loc[df["library_layout"] == "PAIRED", "read_length"] /= 2
    df["chemistry"] = df["library_selection"].map(cfg["library_selection"])
    df["canonical_tissue"] = df["tissue"].map(lambda t: cfg["tissues"][t]["canonical"])
    df["mt_content"] = df["tissue"].map(lambda t: cfg["tissues"][t]["mt_content"])
    text = (df["study_title"].fillna("") + " " + df["sample_title"].fillna("")
            + " " + df["sample_description"].fillna("")
            + " " + df["description"].fillna("")).str.lower()
    # "Non-tumor RNA sample" must not trip the tumour filter, so require that the
    # term is not immediately preceded by a negation.
    pat = "|".join(t.lower() for t in cfg["exclude_title_terms"])
    df["excluded_term"] = text.str.extract(f"(?<!non-)(?<!non )({pat})", expand=False)
    ppat = "|".join(t.lower() for t in cfg["exclude_protocol_terms"])
    df["excluded_protocol"] = text.str.extract(f"({ppat})", expand=False)
    return df


def summarise(df, cfg):
    """Collapse runs to one row per (study, tissue), then to study candidates."""
    keep = df[df["excluded_term"].isna() & df["excluded_protocol"].isna()
              & (df["read_length"] >= cfg["min_read_length"])]
    # A run can match several keywords (brain AND cortex); count it once.
    keep = keep.drop_duplicates(["run_accession", "canonical_tissue"])
    per_tissue = (
        keep.groupby(["study_accession", "canonical_tissue"])
        .agg(n_runs=("run_accession", "count"),
             n_samples=("sample_accession", "nunique"),
             median_read_length=("read_length", "median"),
             median_read_count=("read_count", "median"),
             layouts=("library_layout", lambda s: "/".join(sorted(set(s.dropna())))),
             chemistries=("chemistry", lambda s: "/".join(sorted(set(s.dropna())))),
             mt_content=("mt_content", "first"),
             study_title=("study_title", "first"))
        .reset_index()
    )
    per_tissue = per_tissue[per_tissue["n_runs"] >= cfg["min_runs_per_tissue"]]

    out = (
        per_tissue.groupby("study_accession")
        .agg(n_tissues=("canonical_tissue", "nunique"),
             tissues=("canonical_tissue", lambda s: ",".join(sorted(set(s)))),
             mt_content_span=("mt_content", lambda s: ",".join(sorted(set(s.dropna())))),
             n_runs=("n_runs", "sum"),
             n_samples=("n_samples", "sum"),
             median_read_length=("median_read_length", "median"),
             median_read_count=("median_read_count", "median"),
             layouts=("layouts", lambda s: "/".join(sorted(set("/".join(s).split("/"))))),
             chemistries=("chemistries", lambda s: "/".join(sorted(set("/".join(s).split("/"))))),
             study_title=("study_title", "first"))
        .reset_index()
    )
    out = out[(out["n_tissues"] >= cfg["min_tissues_per_study"])
              & (out["n_runs"] >= cfg["min_runs_per_study"])]
    # Rank by what the design needs: tissue span, then chemistry purity, then depth.
    out["role_hint"] = out["chemistries"].map(
        {"polya": "primary", "rrna_depleted": "replication"}).fillna("mixed/unclear")
    out["span_score"] = out["mt_content_span"].str.count(",") + 1
    return out.sort_values(["span_score", "n_runs"], ascending=False)


def selftest():
    cfg = {"library_selection": {"PolyA": "polya"}, "tissues": {"heart": {"canonical": "heart", "mt_content": "high"},
                       "blood": {"canonical": "blood", "mt_content": "low"}},
           "exclude_title_terms": ["tumor"], "exclude_protocol_terms": ["snrna"],
           "min_read_length": 50,
           "min_runs_per_tissue": 2, "min_tissues_per_study": 2, "min_runs_per_study": 3}
    rows = []
    for i in range(3):  # heart, healthy, paired 2x100
        rows.append(dict(study_accession="P1", sample_accession=f"s{i}", run_accession=f"r{i}",
                         instrument_model="X", library_layout="PAIRED", library_selection="PolyA",
                         read_count="1000", base_count="200000", study_title="healthy heart",
                         sample_title="heart", sample_description="", tissue_type="heart",
                         description="", first_public="", tissue="heart"))
    for i in range(2):
        rows.append(dict(study_accession="P1", sample_accession=f"b{i}", run_accession=f"q{i}",
                         instrument_model="X", library_layout="PAIRED", library_selection="PolyA",
                         read_count="1000", base_count="200000", study_title="healthy heart",
                         sample_title="blood", sample_description="", tissue_type="blood",
                         description="", first_public="", tissue="blood"))
    rows.append(dict(study_accession="P2", sample_accession="t0", run_accession="t0",
                     instrument_model="X", library_layout="PAIRED", library_selection="PolyA",
                     read_count="1000", base_count="200000", study_title="tumor study",
                     sample_title="heart", sample_description="", tissue_type="heart",
                         description="", first_public="", tissue="heart"))
    a = annotate(pd.DataFrame(rows), cfg)
    assert a["read_length"].iloc[0] == 100, a["read_length"].iloc[0]   # paired halving
    assert a["excluded_term"].iloc[-1] == "tumor"                      # healthy filter bites
    neg = annotate(pd.DataFrame([dict(rows[0], study_title="Non-tumor heart sample",
                                      sample_title="heart", description="snRNA-seq")]), cfg)
    assert pd.isna(neg["excluded_term"].iloc[0])       # "non-tumor" is not a tumour study
    assert neg["excluded_protocol"].iloc[0] == "snrna" # but single-cell is still dropped
    s = summarise(a, cfg)
    assert list(s["study_accession"]) == ["P1"], list(s["study_accession"])
    assert s["n_tissues"].iloc[0] == 2 and s["span_score"].iloc[0] == 2
    assert s["role_hint"].iloc[0] == "primary"
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/params.yaml")
    ap.add_argument("--out", default="results/tables/bioproject_candidates.tsv")
    ap.add_argument("--runs-out", default="results/tables/bioproject_scan_runs.tsv.gz")
    ap.add_argument("--raw-out", default="results/tables/bioproject_scan_raw.tsv.gz")
    ap.add_argument("--from-raw", action="store_true",
                    help="re-summarise the cached raw pull instead of hitting ENA")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        return selftest()

    cfg = _config.load(args.config)["bioproject_scan"]
    if args.from_raw:
        runs = annotate(pd.read_csv(args.raw_out, sep="\t", dtype=str), cfg)
        cand = summarise(runs, cfg)
        cand.to_csv(args.out, sep="\t", index=False)
        print(f"{len(cand)} candidate studies -> {args.out}")
        return
    frames = []
    for selection in cfg["library_selection"]:
        for tissue in cfg["tissues"]:
            df = ena_query(cfg, selection, tissue)
            print(f"{selection:20s} {tissue:16s} {len(df):7d} runs", file=sys.stderr)
            frames.append(df)
    runs = pd.concat(frames, ignore_index=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    runs.to_csv(args.raw_out, sep="\t", index=False)
    runs = annotate(runs, cfg)
    runs.to_csv(args.runs_out, sep="\t", index=False)
    cand = summarise(runs, cfg)
    cand.to_csv(args.out, sep="\t", index=False)
    print(f"\n{len(cand)} candidate studies -> {args.out}")
    print("STOP. Cohort selection is a human decision (AGENTS.md 7.2).")


if __name__ == "__main__":
    main()
