#!/usr/bin/env python3
"""Phase 6b - compare RNA-derived haplogroups against the donors' own DNA.

GEUVADIS donors are 1000 Genomes individuals, so the same people have a
whole-genome chrM callset. Both sides are classified with the same Haplogrep
version from the same reference, so a disagreement is the RNA pipeline's, not
the tool's.

Also answers the question Phase 4 raised: the control region is the shallowest
part of the molecule in RNA-seq, and it is where haplogroup-defining variants
concentrate - so are the positions that actually drive assignment covered?
"""
import argparse, gzip, subprocess, sys, urllib.request
from pathlib import Path
import pandas as pd
import _config


def fetch_truth_vcf(url, dest):
    dest = Path(dest)
    if not dest.exists():
        print(f"fetching {dest.name}", file=sys.stderr)
        urllib.request.urlretrieve(url, dest)
    return dest


def subset_vcf(vcf, wanted, out):
    """Keep only our donors' columns. Written by hand because the container has
    no bcftools, and pulling one in for a column subset is not worth a rebuild."""
    keep_idx, header = None, []
    with gzip.open(vcf, "rt") as f, open(out, "w") as g:
        for line in f:
            if line.startswith("##"):
                g.write(line)
                continue
            p = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                keep_idx = [i for i, s in enumerate(p[9:], start=9) if s in wanted]
                found = [p[i] for i in keep_idx]
                missing = sorted(wanted - set(found))
                if missing:
                    print(f"not in the 1000G callset: {', '.join(missing)}", file=sys.stderr)
                header = p[:9] + found
                g.write("\t".join(header) + "\n")
                continue
            g.write("\t".join(p[:9] + [p[i] for i in keep_idx]) + "\n")
    return len(header) - 9


def classify(path, out, fmt="vcf"):
    subprocess.run(["haplogrep", "classify", "--format", fmt, "--extend-report",
                    "--in", str(path), "--out", str(out)], check=True)
    return pd.read_csv(out, sep="\t")


def trim(hg, depth):
    """Haplogroup truncated to `depth` tree levels, for graded agreement.

    Exact string equality is the wrong bar on its own: H1b1 and H1b differ by one
    subclade split, which a missing low-coverage position can cause, and calling
    that a failure would hide how close the pipeline actually is.
    """
    if not isinstance(hg, str) or not hg:
        return ""
    out, letters = "", 0
    for ch in hg:
        out += ch
        if ch.isdigit() or ch.isalpha():
            letters += 1
        if letters >= depth:
            break
    return out[:depth]


def compare(rna, dna, samples):
    """RNA vs DNA per donor, and A vs B within RNA."""
    rna = rna.copy()
    rna["run_accession"] = rna["SampleID"].str.rsplit("_", n=1).str[0]
    rna["arm"] = rna["SampleID"].str.rsplit("_", n=1).str[1]
    m = samples[["run_accession", "sample_title", "bioproject"]]
    rna = rna.merge(m, on="run_accession", how="left")

    arms = rna.pivot_table(index="run_accession", columns="arm", values="Haplogroup",
                           aggfunc="first").reset_index()
    arms["A_equals_B"] = arms["A"] == arms["B"]

    d = dna.rename(columns={"SampleID": "sample_title", "Haplogroup": "hg_dna",
                            "Quality": "q_dna"})
    j = (rna[rna.arm == "A"]
         .rename(columns={"Haplogroup": "hg_rna", "Quality": "q_rna"})
         .merge(d[["sample_title", "hg_dna", "q_dna"]], on="sample_title", how="inner"))
    j["exact"] = j.hg_rna == j.hg_dna
    for k in (1, 2, 3):
        j[f"match_{k}"] = j.hg_rna.map(lambda x: trim(x, k)) == j.hg_dna.map(lambda x: trim(x, k))
    return arms, j


def selftest():
    assert trim("H1b1", 1) == "H" and trim("H1b1", 2) == "H1" and trim("H1b1", 3) == "H1b"
    assert trim("L3e1a1", 1) == "L"
    assert trim("", 2) == "" and trim(None, 2) == ""
    # A deeper request than the label has must return the whole label, not pad it.
    assert trim("H", 3) == "H"

    rna = pd.DataFrame({"SampleID": ["r1_A", "r1_B", "r2_A", "r2_B"],
                        "Haplogroup": ["H1b1", "H1b1", "U5a", "U5a1"],
                        "Quality": [.95, .95, .9, .9]})
    dna = pd.DataFrame({"SampleID": ["NA1", "NA2"], "Haplogroup": ["H1b1", "U5b"],
                        "Quality": [1.0, 1.0]})
    s = pd.DataFrame({"run_accession": ["r1", "r2"], "sample_title": ["NA1", "NA2"],
                      "bioproject": ["P", "P"]})
    arms, j = compare(rna, dna, s)
    assert arms.set_index("run_accession").A_equals_B.tolist() == [True, False]
    assert j.set_index("run_accession").exact.tolist() == [True, False]
    # U5a vs U5b agree at two levels and part at the third.
    r2 = j.set_index("run_accession").loc["r2"]
    assert r2.match_1 and r2.match_2 and not r2.match_3, r2.to_dict()
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
    h = c["haplogroup"]
    hdir, out = Path(h["hsd_dir"]), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(a.samples, sep="\t")

    rna = pd.read_csv(hdir / "haplogroups_rna.txt", sep="\t")
    vcf = fetch_truth_vcf(h["truth_vcf_url"], hdir / "1000G_chrMT.vcf.gz")
    wanted = set(samples[samples.bioproject == "PRJEB3366"].sample_title.astype(str))
    n = subset_vcf(vcf, wanted, hdir / "truth_subset.vcf")
    print(f"{n} donors found in the 1000G chrM callset", file=sys.stderr)
    dna = classify(hdir / "truth_subset.vcf", hdir / "haplogroups_dna.txt")

    arms, j = compare(rna, dna, samples)
    arms.to_csv(out / "haplogroup_arm_agreement.tsv", sep="\t", index=False)
    j.to_csv(out / "haplogroup_rna_vs_dna.tsv", sep="\t", index=False)
    print(f"A==B in {arms.A_equals_B.sum()}/{len(arms)} samples")
    print(f"RNA==DNA exactly in {j.exact.sum()}/{len(j)} donors; "
          f"first level {j.match_1.sum()}/{len(j)}, "
          f"two levels {j.match_2.sum()}/{len(j)}")


if __name__ == "__main__":
    main()
