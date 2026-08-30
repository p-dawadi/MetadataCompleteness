#!/usr/bin/env python3
"""
type_strain_identity_check.py

A real, alignment-based substitute for BLAST/type-strain verification
(Reviewer 3, comment R3-8: "Species identification was not independently
confirmed"). NCBI/BLAST web access was not available in this environment
(robots.txt/CAPTCHA blocked), so instead of leaving this open we used data
already in hand: the curated dataset happens to include the *S. anginosus*
type strain itself, twice over --

  - AF104678.1  "Streptococcus anginosus strain ATCC33397 16S ribosomal RNA
                 gene, partial sequence"   (ATCC 33397 = the type strain)
  - NR_041722.2 "Streptococcus anginosus SK52 = DSM 20563 16S ribosomal RNA,
                 complete sequence"        (an NCBI RefSeq targeted-locus
                 (NR_) record -- NCBI's own curated representative/type
                 material sequence; DSM 20563 is the same type strain
                 deposited under a different culture-collection number)

Both are already part of the final alignment (Aligned.Strepto.aln). This
script computes, directly from that alignment (no external network access
required), the ungapped pairwise percent identity of every other curated
sequence against each type-strain reference, and reports the lowest-identity
records as the ones that most need independent BLAST/phylogenetic-placement
follow-up.

Usage:
    python type_strain_identity_check.py [--threshold 95.0]

Requires: none beyond the standard library (reads FASTA directly).
"""
import argparse
import re

ALIGNMENT = "../alignment/Aligned.Strepto.aln"
REFERENCES = {
    "AF104678.1": "ATCC 33397 (type strain)",
    "NR_041722.2": "SK52 = DSM 20563 (RefSeq targeted-locus type-strain record)",
}


def read_fasta(path):
    seqs = {}
    label = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if label is not None:
                    seqs[label] = "".join(chunks)
                label = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if label is not None:
        seqs[label] = "".join(chunks)
    return seqs


def pct_identity(a, b):
    """Ungapped percent identity over aligned columns where neither
    sequence has a gap (standard pairwise-identity-from-alignment
    definition)."""
    compared = 0
    matches = 0
    for x, y in zip(a, b):
        if x in "-.Nn" or y in "-.Nn":
            continue
        compared += 1
        if x.upper() == y.upper():
            matches += 1
    if compared == 0:
        return None, 0
    return 100.0 * matches / compared, compared


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threshold", type=float, default=95.0,
                     help="Flag records below this %% identity to the type strain (default 95.0).")
    args = ap.parse_args()

    seqs = read_fasta(ALIGNMENT)
    print(f"Alignment: {len(seqs)} sequences, {len(next(iter(seqs.values())))} columns\n")

    for ref_label, ref_desc in REFERENCES.items():
        if ref_label not in seqs:
            print(f"[skip] {ref_label} not found in alignment")
            continue
        ref_seq = seqs[ref_label]
        print(f"=== Reference: {ref_label} ({ref_desc}) ===")
        results = []
        for label, seq in seqs.items():
            if label == ref_label:
                continue
            ident, n = pct_identity(ref_seq, seq)
            if ident is not None:
                results.append((ident, n, label))
        results.sort()

        n_below = sum(1 for ident, n, lbl in results if ident < args.threshold)
        idents = [r[0] for r in results]
        print(f"n compared: {len(results)}  mean identity: {sum(idents)/len(idents):.1f}%  "
              f"median: {sorted(idents)[len(idents)//2]:.1f}%  "
              f"min: {min(idents):.1f}%  max: {max(idents):.1f}%")
        print(f"records below {args.threshold:.1f}% identity to this reference: {n_below}")
        print("10 lowest-identity records (accession, %identity, compared columns):")
        for ident, n, lbl in results[:10]:
            print(f"  {lbl:20s} {ident:6.1f}%  (n={n})")
        print()


if __name__ == "__main__":
    main()
