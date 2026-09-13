#!/usr/bin/env python3
"""Phase 6 - haplogroup assignment from RNA-seq, checked against DNA.

Three things, in order:

  call      pileup one sample's chrM BAM and write a Haplogrep .hsd profile
  classify  run Haplogrep over the collected profiles
  compare   RNA vs DNA for the GEUVADIS donors, A vs B, and the coverage of the
            positions Haplogrep actually used

The variant call here is a *majority* call. Haplogroups are defined by near-fixed
variants, so a low-frequency threshold would let an RT artefact or a NUMT read
rename a lineage. Low-frequency detection is Phase 7's subject, not this one.
"""
import argparse, subprocess, sys
from collections import Counter
from pathlib import Path
import pandas as pd
import yaml

BASES = "ACGT"


def read_reference(path):
    seq = []
    for line in open(path):
        if not line.startswith(">"):
            seq.append(line.strip())
    return "".join(seq).upper()


def call_variants(counts, ref, cfg):
    """Pileup base counts -> Haplogrep polymorphism strings.

    `counts` maps 1-based position to a Counter of observed bases. A position is
    called only when it is deep enough and one non-reference base holds at least
    af_threshold of the reads; everything else is left uncalled rather than
    guessed, so Haplogrep sees a gap instead of a wrong base.
    """
    skip = set(cfg["skip_positions"])
    polys, called, undercovered = [], 0, 0
    for pos in range(1, len(ref) + 1):
        if pos in skip:
            continue
        c = counts.get(pos)
        depth = sum(c[b] for b in BASES) if c else 0
        if depth < cfg["min_depth"]:
            undercovered += 1
            continue
        called += 1
        top, n = max(((b, c[b]) for b in BASES), key=lambda x: x[1])
        if top != ref[pos - 1] and n / depth >= cfg["af_threshold"]:
            polys.append(f"{pos}{top}")
    return polys, called, undercovered


def pileup_counts(bam, length, min_bq):
    """Per-position base counts from a BAM, ignoring secondary/supplementary."""
    import pysam
    counts = {}
    with pysam.AlignmentFile(bam, "rb") as b:
        for col in b.pileup(max_depth=1_000_000, min_base_quality=min_bq,
                            ignore_overlaps=False, stepper="samtools"):
            pos = col.reference_pos + 1
            if pos > length:
                continue
            c = Counter()
            for base in col.get_query_sequences(add_indels=False):
                base = base.upper()
                if base in BASES:
                    c[base] += 1
            if c:
                counts[pos] = c
    return counts


def hsd_line(sample, length, polys):
    return f"{sample}\t1-{length}\t?\t{' '.join(polys)}\n"


def do_call(row, cfg, ref_path, bam_dir, out_dir):
    ref = read_reference(ref_path)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    run = row["run_accession"]
    rows = []
    for arm in ("A", "B"):
        bam = Path(bam_dir) / f"{run}_{arm}_chrM.bam"
        counts = pileup_counts(bam, len(ref), cfg["min_base_quality"])
        polys, called, under = call_variants(counts, ref, cfg)
        (out_dir / f"{run}_{arm}.hsd").write_text(
            "SampleId\tRange\tHaplogroup\tPolymorphisms\n"
            + hsd_line(f"{run}_{arm}", len(ref), polys))
        rows.append({"run_accession": run, "arm": arm, "n_polys": len(polys),
                     "positions_called": called, "positions_undercovered": under})
        print(f"{run} {arm}: {len(polys)} polymorphisms, "
              f"{under} positions below {cfg['min_depth']}x", file=sys.stderr)
    pd.DataFrame(rows).to_csv(out_dir / f"{run}_callstats.tsv", sep="\t", index=False)


def do_classify(out_dir, sif=None):
    """Concatenate the per-sample profiles and classify them in one Haplogrep run."""
    out_dir = Path(out_dir)
    files = sorted(out_dir.glob("*_[AB].hsd"))
    if not files:
        sys.exit(f"no .hsd profiles in {out_dir}")
    combined = out_dir / "all_samples.hsd"
    with open(combined, "w") as f:
        f.write("SampleId\tRange\tHaplogroup\tPolymorphisms\n")
        for p in files:
            f.write(p.read_text().splitlines()[1] + "\n")
    out = out_dir / "haplogroups_rna.txt"
    cmd = ["haplogrep", "classify", "--format", "hsd", "--extend-report",
           "--in", str(combined), "--out", str(out)]
    subprocess.run(cmd, check=True)
    print(f"classified {len(files)} profiles -> {out}")


def selftest():
    ref = "ACGTACGTAC"
    cfg = {"min_depth": 10, "af_threshold": 0.9, "skip_positions": [5]}
    counts = {
        1: Counter({"A": 20}),                     # matches reference
        2: Counter({"T": 19, "C": 1}),             # 95% alt -> called
        3: Counter({"T": 12, "G": 8}),             # 60% alt -> not fixed, skipped
        4: Counter({"A": 5}),                      # too shallow
        5: Counter({"G": 40}),                     # placeholder position, skipped
        6: Counter({"C": 30}),                     # matches reference (ref[5] == C)
    }
    polys, called, under = call_variants(counts, ref, cfg)
    assert polys == ["2T"], polys
    assert called == 4, called                     # 1,2,3,6 -- not 4 (shallow), not 5
    assert under == 5, under                       # 4,7,8,9,10
    # A skipped position must not be silently called even at huge depth.
    assert not any(p.startswith("5") for p in polys)

    line = hsd_line("s1", 16569, ["263G", "16519C"])
    assert line == "s1\t1-16569\t?\t263G 16519C\n", repr(line)
    assert hsd_line("s1", 16569, []).endswith("?\t\n"), "empty profile must still be a row"

    import tempfile, os
    f = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
    f.write(">chrM\nACGT\nACGT\n"); f.close()
    assert read_reference(f.name) == "ACGTACGT"
    os.unlink(f.name)
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mode", choices=["call", "classify", "selftest"])
    ap.add_argument("--config", default="config/params.yaml")
    ap.add_argument("--samples", default="config/samples.tsv")
    ap.add_argument("--row", type=int)
    a = ap.parse_args()
    if a.mode == "selftest":
        return selftest()
    c = yaml.safe_load(open(a.config))
    hcfg, acfg = c["haplogroup"], c["alignment"]
    if a.mode == "call":
        if not a.row:
            sys.exit("--row is required for call")
        row = pd.read_csv(a.samples, sep="\t", dtype=str).iloc[a.row - 1]
        do_call(row, hcfg, Path(acfg["ref_dir"]) / "chrM_rCRS.fa",
                acfg["bam_dir"], hcfg["hsd_dir"])
    else:
        do_classify(hcfg["hsd_dir"])


if __name__ == "__main__":
    main()
