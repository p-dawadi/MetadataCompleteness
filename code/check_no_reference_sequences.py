#!/usr/bin/env python3
"""
check_no_reference_sequences.py

Confirms that the final alignment/tree (Aligned.Strepto.aln, 354 sequences)
contains no external reference sequences of any kind -- i.e., no SILVA
database sequences themselves are present among the 354 aligned/tree tips.
This was originally run to check for leftover SILVA reference sequences
after the manuscript's original submission mentioned SILVA; the authors
subsequently confirmed that SILVA genuinely was used, as a screening step to
confirm 16S rRNA gene identity among candidate sequences prior to the final
MAFFT alignment (not as the MAFFT alignment reference itself, and not by
including SILVA's own reference sequences in the final aligned set). This
script's result is unaffected by that clarification and remains a useful,
independent check: every tip label matches either the 350 curated study
accessions (data/curated_metadata.csv) or the four additional GenBank S.
anginosus records identified independently during tree-metadata
cross-referencing (PX419550, PX419553, PX419555, PX419685). Any label that
matches neither would indicate an external or unexplained sequence; this
script reports that count (expected: zero).

Usage: python check_no_reference_sequences.py
"""
import re
import pandas as pd

ALIGNMENT = "../alignment/Aligned.Strepto.aln"
METADATA = "../data/curated_metadata.csv"
KNOWN_EXTRA = {"PX419550", "PX419553", "PX419555", "PX419685"}


def base_accession(label):
    return re.sub(r"\.\d+$", "", label)


def main():
    df = pd.read_csv(METADATA)
    curated_base = set(base_accession(a.strip()) for a in df["Accession Number"].astype(str))

    labels = []
    with open(ALIGNMENT) as fh:
        for line in fh:
            if line.startswith(">"):
                labels.append(line[1:].split()[0])

    base_labels = [base_accession(l) for l in labels]
    unexplained = [l for l, b in zip(labels, base_labels) if b not in curated_base and b not in KNOWN_EXTRA]

    print(f"Alignment/tree tip labels: {len(labels)}")
    print(f"Curated study accessions: {len(curated_base)}")
    print(f"Matched to curated dataset: {sum(1 for b in base_labels if b in curated_base)}")
    print(f"Matched to known additional GenBank records: {sum(1 for b in base_labels if b in KNOWN_EXTRA)}")
    print(f"Unexplained labels (would indicate leftover reference/other sequences): {len(unexplained)}")
    for u in unexplained:
        print(f"  {u}")
    if not unexplained:
        print("\nResult: every tip is accounted for. No external reference "
              "sequences (SILVA or otherwise) are themselves present as tips "
              "in the final alignment/tree -- consistent with an alignment "
              "built from the retrieved sequences, with SILVA used only as a "
              "screening step to confirm 16S rRNA gene identity beforehand "
              "(see manuscript Methods: Sequence alignment), not as a source "
              "of added reference sequences.")


if __name__ == "__main__":
    main()
