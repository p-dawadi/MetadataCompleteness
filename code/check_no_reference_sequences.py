#!/usr/bin/env python3
"""
check_no_reference_sequences.py

Confirms that every tip label in an alignment/tree traces back to a real
curated study accession (data/curated_metadata.csv) and that no external
reference or database sequence (SILVA or otherwise) has been left in the
aligned/tree set by accident. Any tip label that does not match a curated
accession is reported explicitly (expected count: zero, now that the 4
additional non-curated GenBank records, PX419550/553/555/685, have been
removed from the dataset entirely).

Usage:
    python check_no_reference_sequences.py Strepto_350_curated_only.fasta
    python check_no_reference_sequences.py Aligned_346_sensitivity.fasta
"""
import argparse
import re
import pandas as pd


def base_accession(label):
    return re.sub(r"\.\d+$", "", label)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", help="Alignment or unaligned FASTA whose tip/sequence labels to check")
    ap.add_argument("--metadata", default="curated_metadata.csv",
                     help="Path to curated_metadata.csv (default: curated_metadata.csv in the current dir)")
    args = ap.parse_args()

    df = pd.read_csv(args.metadata)
    curated_base = set(base_accession(a.strip()) for a in df["Accession Number"].astype(str))

    labels = []
    with open(args.fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                labels.append(line[1:].split()[0])

    base_labels = [base_accession(l) for l in labels]
    unexplained = [l for l, b in zip(labels, base_labels) if b not in curated_base]

    print(f"Alignment/tree tip labels: {len(labels)}")
    print(f"Curated study accessions: {len(curated_base)}")
    print(f"Matched to curated dataset: {sum(1 for b in base_labels if b in curated_base)}")
    print(f"Unexplained labels (would indicate leftover reference/other sequences): {len(unexplained)}")
    for u in unexplained:
        print(f"  {u}")
    if not unexplained:
        print("\nResult: every tip is accounted for against the curated 350-record "
              "dataset. No external reference sequences (SILVA or otherwise) are "
              "themselves present as tips, and no non-curated GenBank record "
              "remains in the set.")


if __name__ == "__main__":
    main()
