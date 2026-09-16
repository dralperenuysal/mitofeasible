#!/usr/bin/env python3
"""Phase 3a - build every reference the alignment needs.

Produces, under ref/:
  genome_A.fa        GRCh38 primary assembly, NUMTs left alone
  genome_B.fa        same, NUMTs hard-masked to N in the nuclear genome only
  numt_regions.bed   UCSC hg38 nuMtSeq track
  chrM_rCRS.fa       chrM alone
  chrM_shifted.fa    chrM rotated by `reference.shift`, for the control region
  annotation.gtf     GENCODE, for STAR splice junctions

A and B differ in the reference and nothing else (AGENTS.md 3, Phase 3).
"""
import argparse, gzip, shutil, sys, urllib.request
from pathlib import Path
import _config

GENCODE = ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/"
           "GRCh38.primary_assembly.genome.fa.gz")
GTF = ("https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/"
       "gencode.v47.primary_assembly.annotation.gtf.gz")
NUMT_API = ("https://api.genome.ucsc.edu/getData/track?genome=hg38;track=nuMtSeq;"
            "maxItemsOutput=100000")


def fetch(url, dest):
    if dest.exists():
        print(f"have {dest.name}", file=sys.stderr)
        return dest
    print(f"fetching {dest.name}", file=sys.stderr)
    tmp = dest.with_suffix(dest.suffix + ".part")
    urllib.request.urlretrieve(url, tmp)
    tmp.rename(dest)
    return dest


def gunzip(src, dest):
    if dest.exists():
        return dest
    with gzip.open(src, "rb") as f, open(dest, "wb") as g:
        shutil.copyfileobj(f, g)
    return dest


def numt_bed(dest):
    """UCSC hg38 nuMtSeq track as BED3. Fetched, not vendored, so the provenance
    is a URL rather than a file of unknown age."""
    if dest.exists():
        return dest
    import json
    with urllib.request.urlopen(NUMT_API, timeout=300) as r:
        rows = json.load(r)["nuMtSeq"]
    if isinstance(rows, dict):            # UCSC returns per-chrom dicts when large
        rows = [x for v in rows.values() for x in v]
    with open(dest, "w") as f:
        for x in sorted(rows, key=lambda r: (r["chrom"], r["chromStart"])):
            f.write(f'{x["chrom"]}\t{x["chromStart"]}\t{x["chromEnd"]}\t{x["name"]}\n')
    print(f"{len(rows)} NUMT intervals", file=sys.stderr)
    return dest


def read_fasta(path):
    """Yield (header, sequence) pairs. The genome is 3 GB, so one record at a time."""
    name, chunks = None, []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def write_fasta(handle, name, seq, width=60):
    handle.write(f">{name}\n")
    for i in range(0, len(seq), width):
        handle.write(seq[i:i + width] + "\n")


def mask(seq, intervals):
    """Hard-mask half-open BED intervals to N."""
    s = list(seq)
    for start, end in intervals:
        for i in range(max(0, start), min(len(s), end)):
            s[i] = "N"
    return "".join(s)


def rotate(seq, shift):
    """Rotate a circular sequence left by `shift`: new position 1 is old shift+1."""
    shift %= len(seq)
    return seq[shift:] + seq[:shift]


def build(ref_dir, cfg):
    ref_dir.mkdir(parents=True, exist_ok=True)
    fa_gz = fetch(GENCODE, ref_dir / "GRCh38.primary_assembly.genome.fa.gz")
    fetch(GTF, ref_dir / "annotation.gtf.gz")
    gunzip(ref_dir / "annotation.gtf.gz", ref_dir / "annotation.gtf")
    bed = numt_bed(ref_dir / "numt_regions.bed")
    fa = gunzip(fa_gz, ref_dir / "genome_A.fa")

    by_chrom = {}
    for line in open(bed):
        c, s, e = line.split("\t")[:3]
        by_chrom.setdefault(c, []).append((int(s), int(e)))

    chrM = None
    with open(ref_dir / "genome_B.fa", "w") as out:
        for header, seq in read_fasta(fa):
            name = header.split()[0]
            if name == "chrM":
                chrM = seq
                # Never mask chrM itself: the NUMT track marks nuclear copies of
                # mtDNA, so masking chrM would delete the study's subject.
                write_fasta(out, header, seq)
            else:
                write_fasta(out, header, mask(seq, by_chrom.get(name, [])))
    if chrM is None:
        sys.exit("no chrM in the assembly")
    if len(chrM) != cfg["length"]:
        sys.exit(f"chrM is {len(chrM)} bp, expected {cfg['length']} (not rCRS?)")

    with open(ref_dir / "chrM_rCRS.fa", "w") as f:
        write_fasta(f, "chrM", chrM)
    with open(ref_dir / "chrM_shifted.fa", "w") as f:
        write_fasta(f, f"chrM_shifted_{cfg['shift']}", rotate(chrM, cfg["shift"]))
    print(f"ok: chrM {len(chrM)} bp, shifted by {cfg['shift']}")


def selftest():
    seq = "AAAACCCCGGGGTTTT"
    assert mask(seq, [(4, 8)]) == "AAAANNNNGGGGTTTT"
    assert mask(seq, [(-5, 2), (14, 99)]) == "NNAACCCCGGGGTTNN", "clamping"
    assert mask(seq, []) == seq
    r = rotate(seq, 4)
    assert r == "CCCCGGGGTTTTAAAA", r
    # Rotating back must restore the original: the coordinate fold in Phase 4
    # depends on this being exact.
    assert rotate(r, len(seq) - 4) == seq
    assert rotate(seq, len(seq)) == seq
    import io
    h = io.StringIO()
    write_fasta(h, "x", "ACGT" * 20, width=60)
    body = "".join(h.getvalue().split("\n")[1:])
    assert body == "ACGT" * 20, "write_fasta lost bases"
    print("selftest ok")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config/params.yaml")
    ap.add_argument("--ref-dir", default="ref")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        return selftest()
    build(Path(a.ref_dir), _config.load(a.config)["reference"])


if __name__ == "__main__":
    main()
