#!/usr/bin/env python3
"""Phase 2 - build the sample manifest and fetch raw FASTQ for the two cohorts.

Two steps, deliberately separate: --manifest writes config/samples.tsv from ENA
metadata and stops; --download fetches the FASTQ named there and verifies the
ENA md5 of every file. Raw FASTQ only - no count matrices (AGENTS.md 3, Phase 2).
"""
import argparse, hashlib, pathlib, subprocess, sys, urllib.parse, urllib.request
from pathlib import Path
import pandas as pd
import yaml

ENA = "https://www.ebi.ac.uk/ena/portal/api/search"
FIELDS = [
    "run_accession", "sample_accession", "sample_alias", "sample_title",
    "library_selection", "library_layout", "read_count", "base_count",
    "instrument_model", "fastq_ftp", "fastq_md5", "fastq_bytes",
]


def ena_runs(study):
    url = ENA + "?" + urllib.parse.urlencode({
        "result": "read_run", "query": f'study_accession="{study}"',
        "fields": ",".join(FIELDS), "format": "tsv", "limit": 10000,
    })
    with urllib.request.urlopen(url, timeout=300) as r:
        return pd.read_csv(pd.io.common.StringIO(r.read().decode()), sep="\t", dtype=str)


def panel(cfg):
    """1000 Genomes sample panel: the only source of sex/population for GEUVADIS."""
    p = Path(cfg["panel_cache"])
    if not p.exists():
        p.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(cfg["panel_url"], p)
    return pd.read_csv(p, sep="\t", dtype=str).dropna(axis=1, how="all")


def pick_runs(df, spec, seed, pan):
    """One run per sample, then a seeded balanced subsample of n_target samples."""
    df = df.copy()
    df["read_count"] = pd.to_numeric(df["read_count"], errors="coerce")
    df["base_count"] = pd.to_numeric(df["base_count"], errors="coerce")
    df["read_length"] = (df["base_count"] / df["read_count"]).round()
    df.loc[df["library_layout"] == "PAIRED", "read_length"] /= 2
    if spec["panel_join"]:
        # Inner join on purpose: a GEUVADIS LCL absent from the phase 3 panel has no
        # WGS truth set, and matched DNA is the reason this cohort was chosen.
        df = df.merge(pan, left_on="sample_title", right_on="sample", how="inner")
        df["sex"] = df["gender"]
        df["age_class"] = "adult"
    elif spec.get("meta_tsv"):
        # Age, sex and timepoint are in the GEO record, not in any ENA field;
        # the curated join lives in ref/ so the selection stays reproducible.
        meta = pd.read_csv(spec["meta_tsv"], sep="\t", dtype=str)
        df = df.merge(meta, on="run_accession", how="inner")
        for col, val in (spec.get("keep_where") or {}).items():
            df = df[df[col] == val]
    else:
        df["sex"] = pd.NA
        df["age_class"] = pd.NA

    # Deduplicate after filtering, not before: otherwise a run that the metadata
    # filter is about to drop can win the dedup and take its whole sample with it.
    # Technical replicates exist (GEUVADIS runs the same LCL in two labs); keep the
    # deepest run rather than merging, which would mix batch into depth.
    df = df.sort_values("read_count", ascending=False).drop_duplicates("sample_accession")

    strata = [c for c in spec["strata"] if c in df.columns]
    df = df.sample(frac=1, random_state=seed)
    # Round-robin across strata so the subsample is balanced without quota maths.
    df["_rr"] = df.groupby(strata, dropna=False).cumcount() if strata else range(len(df))
    df = df.sort_values(["_rr"]).head(spec["n_target"]).drop(columns="_rr")
    return df


def build_manifest(cfg):
    pan = panel(cfg)
    out = []
    for study, spec in cfg["studies"].items():
        d = pick_runs(ena_runs(study), spec, cfg["seed"], pan)
        d["bioproject"] = study
        for k in ("role", "arm", "tissue", "library_type"):
            d[k] = spec[k]
        print(f"{study}: {len(d)} samples selected of {spec['n_target']} requested",
              file=sys.stderr)
        out.append(d)
    cols = ["run_accession", "sample_accession", "sample_title", "bioproject", "role",
            "arm", "tissue", "library_type", "sex", "age_class", "subject_biopsy",
            "super_pop", "pop",
            "library_layout", "read_length", "read_count", "instrument_model",
            "fastq_ftp", "fastq_md5", "fastq_bytes"]
    m = pd.concat(out, ignore_index=True)
    return m[[c for c in cols if c in m.columns]]


def md5(path, chunk=1 << 20):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def fetch_one(url, want, size, f):
    """Download one FASTQ to `f` and verify it. Returns None on success, else why.

    Resume is only safe against a partial file from the same attempt. A file left
    by an aborted earlier run can be truncated mid-record or longer than the
    remote one, and `curl -C -` then either fails ("offset beyond file size") or
    resumes onto garbage and produces a complete-looking file with the wrong md5.
    So: drop anything that is already too big, and on any failure retry once from
    scratch rather than resuming onto a file we no longer trust.
    """
    for attempt, resume in enumerate((True, False)):
        if size and f.exists() and f.stat().st_size > size:
            f.unlink()                       # longer than the source: not a prefix
        if not resume and f.exists():
            f.unlink()
        cmd = ["curl", "-fsSL", "--retry", "10", "--retry-delay", "15",
               "--retry-all-errors", "--speed-time", "60", "--speed-limit", "10000",
               "-o", str(f), "ftp://" + url]
        if resume:
            cmd[2:2] = ["-C", "-"]
        r = subprocess.run(cmd)
        if r.returncode != 0:
            print(f"curl exit {r.returncode} on {f.name}"
                  f"{' (resume)' if resume else ' (fresh)'}", file=sys.stderr)
            continue
        got = md5(f)
        if got == want:
            return None
        print(f"md5 mismatch on {f.name}{' after resume' if resume else ''}",
              file=sys.stderr)
    return f"failed after fresh retry"


def download(manifest, dest):
    """Fetch every FASTQ and verify the ENA md5. Already-verified files are skipped."""
    dest = Path(dest)
    bad = []
    for _, row in manifest.iterrows():
        urls = str(row["fastq_ftp"]).split(";")
        sums = str(row["fastq_md5"]).split(";")
        sizes = [int(x) if x.isdigit() else 0
                 for x in str(row.get("fastq_bytes", "")).split(";")]
        sizes += [0] * (len(urls) - len(sizes))
        d = dest / row["run_accession"]
        d.mkdir(parents=True, exist_ok=True)
        for url, want, size in zip(urls, sums, sizes):
            f = d / url.rsplit("/", 1)[-1]
            if f.exists() and md5(f) == want:
                continue
            print(f"fetching {f.name}", file=sys.stderr)
            why = fetch_one(url, want, size, f)
            if why:
                bad.append((f.name, want, why))
                print(f"FAILED {f.name}: {why}", file=sys.stderr)
    return bad


def selftest():
    pan = pd.DataFrame({"sample": ["NA1", "NA2", "NA3", "NA4"],
                        "pop": ["GBR", "YRI", "GBR", "YRI"],
                        "super_pop": ["EUR", "AFR", "EUR", "AFR"],
                        "gender": ["male", "female", "female", "male"]})
    df = pd.DataFrame({
        "run_accession": list("abcde"),
        "sample_accession": ["S1", "S1", "S2", "S3", "S4"],
        "sample_title": ["NA1", "NA1", "NA2", "NA3", "NA4"],
        "library_layout": ["PAIRED"] * 5,
        "read_count": ["100", "200", "100", "100", "100"],
        "base_count": ["30000", "60000", "30000", "30000", "30000"],
    })
    spec = {"panel_join": True, "strata": ["super_pop", "gender"], "n_target": 2,
            "role": "primary", "arm": "low_mt", "tissue": "lcl", "library_type": "polyA"}
    got = pick_runs(df, spec, 1, pan)
    assert got["super_pop"].notna().all(), "unpanelled sample kept"
    assert got["read_length"].iloc[0] == 150, got["read_length"].tolist()  # paired halved
    assert "a" not in set(got["run_accession"]), "shallower technical replicate kept"
    assert len(got) == 2 and got["super_pop"].nunique() == 2, "strata not balanced"

    import tempfile
    meta = pd.DataFrame({"run_accession": list("abcde"),
                         "age_class": ["young", "young", "old", "old", "young"],
                         "sex": ["male", "male", "female", "male", "female"],
                         "timepoint": ["rest", "3hr_post", "rest", "rest", "rest"]})
    f = tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False)
    meta.to_csv(f.name, sep="\t", index=False)
    spec2 = {"panel_join": False, "strata": ["age_class", "sex"], "n_target": 3,
             "meta_tsv": f.name, "keep_where": {"timepoint": "rest"},
             "role": "replication", "arm": "high_mt", "tissue": "skeletal_muscle",
             "library_type": "polyA"}
    g2 = pick_runs(df, spec2, 1, pan)
    assert "b" not in set(g2["run_accession"]), "non-rest timepoint kept"
    assert "a" in set(g2["run_accession"]), "sample lost because its deeper run was filtered"
    assert g2["age_class"].nunique() == 2, "strata not balanced"
    # fetch_one must delete a local file that is longer than the source before
    # trying to resume onto it - the bug that produced "offset beyond file size".
    import tempfile, os
    d = pathlib.Path(tempfile.mkdtemp())
    f = d / "x.gz"
    f.write_bytes(b"0" * 100)
    calls = []
    real = subprocess.run

    def fake(cmd, *a, **k):
        calls.append(("-C" in cmd, f.exists()))
        f.write_bytes(b"payload")
        return subprocess.CompletedProcess(cmd, 0)

    subprocess.run = fake
    try:
        want = hashlib.md5(b"payload").hexdigest()
        assert fetch_one("h/x.gz", want, 50, f) is None
        assert calls[0] == (True, False), calls   # oversized file removed first
    finally:
        subprocess.run = real

    print("selftest ok")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/params.yaml")
    ap.add_argument("--out", default="config/samples.tsv")
    ap.add_argument("--manifest", action="store_true", help="write samples.tsv and stop")
    ap.add_argument("--download", action="store_true", help="fetch FASTQ from samples.tsv")
    ap.add_argument("--row", type=int, help="1-based row of samples.tsv to fetch "
                    "(for SLURM array jobs; omit to fetch all)")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    cfg = yaml.safe_load(open(a.config))["cohorts"]
    if a.manifest or not a.download:
        m = build_manifest(cfg)
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
        m.to_csv(a.out, sep="\t", index=False)
        gb = pd.to_numeric(m["fastq_bytes"].str.split(";").explode()).sum() / 1e9
        print(f"wrote {a.out}: {len(m)} samples, {gb:.0f} GB to download")
        if not a.download:
            return
    m = pd.read_csv(a.out, sep="\t", dtype=str)
    if a.row:
        m = m.iloc[[a.row - 1]]
    bad = download(m, cfg["fastq_dir"])
    print(f"download complete, {len(bad)} checksum failures")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
