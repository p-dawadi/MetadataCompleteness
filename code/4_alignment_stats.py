#!/usr/bin/env python3
"""
4_alignment_stats.py

Computes the real alignment statistics needed for the manuscript's
"Sequence alignment" Methods paragraph, directly from a MAFFT output
FASTA -- no other tools required (standard library only).

Reports:
  - number of sequences and aligned columns
  - constant / variable / parsimony-informative column counts
  - number of distinct site patterns (columns, deduplicated)
  - gap/missing character % (overall and per-sequence range)
  - ambiguous IUPAC code %
  - mean/median/range of pairwise % identity across all sequence pairs
  - number of unique haplotypes (exact-duplicate aligned sequences collapsed)

Usage:
    python 4_alignment_stats.py Aligned_350_primary.fasta
    python 4_alignment_stats.py Aligned_346_sensitivity.fasta
"""
import argparse
import itertools
import sys
from collections import Counter

GAP_CHARS = set("-.")
AMBIG_CHARS = set("NRYSWKMBDHV")


def read_fasta(path):
    seqs = {}
    label = None
    chunks = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if label is not None:
                    seqs[label] = "".join(chunks)
                label = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
        if label is not None:
            seqs[label] = "".join(chunks)
    return seqs


def pct_identity(a, b):
    compared = matches = 0
    for x, y in zip(a, b):
        if x in GAP_CHARS or y in GAP_CHARS or x in AMBIG_CHARS or y in AMBIG_CHARS:
            continue
        compared += 1
        if x == y:
            matches += 1
    if compared == 0:
        return None
    return 100.0 * matches / compared


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fasta", help="Aligned FASTA file (MAFFT output)")
    ap.add_argument("--max-pairs", type=int, default=None,
                     help="Optional cap on random pairwise comparisons for speed on very large alignments "
                          "(default: all pairs -- fine for a few hundred sequences)")
    ap.add_argument("--seed", type=int, default=20260827)
    args = ap.parse_args()

    seqs = read_fasta(args.fasta)
    labels = list(seqs.keys())
    n_seq = len(labels)
    lengths = {len(s) for s in seqs.values()}
    if len(lengths) != 1:
        sys.exit(f"ERROR: not all sequences are the same length ({lengths}) -- is this really an aligned FASTA?")
    n_col = lengths.pop()
    print(f"Alignment: {n_seq} sequences, {n_col} aligned columns\n")

    # --- column-level stats ---
    constant = variable = pinformative = 0
    site_patterns = Counter()
    total_chars = 0
    gap_chars = 0
    ambig_chars = 0
    per_seq_gap = {lbl: 0 for lbl in labels}

    columns = list(zip(*(seqs[lbl] for lbl in labels)))
    for col in columns:
        pattern = "".join(col)
        site_patterns[pattern] += 1
        non_gap = [c for c in col if c not in GAP_CHARS and c not in AMBIG_CHARS]
        distinct = set(non_gap)
        if len(distinct) <= 1:
            constant += 1
        else:
            variable += 1
            counts = Counter(non_gap)
            n_states_with_ge2 = sum(1 for c in counts.values() if c >= 2)
            if n_states_with_ge2 >= 2:
                pinformative += 1

    for lbl in labels:
        s = seqs[lbl]
        g = sum(1 for c in s if c in GAP_CHARS)
        a = sum(1 for c in s if c in AMBIG_CHARS)
        per_seq_gap[lbl] = 100.0 * g / n_col
        gap_chars += g
        ambig_chars += a
        total_chars += n_col

    gap_pct_overall = 100.0 * gap_chars / total_chars
    ambig_pct_overall = 100.0 * ambig_chars / total_chars
    gap_range = (min(per_seq_gap.values()), max(per_seq_gap.values()))

    print(f"Constant columns: {constant} ({100*constant/n_col:.2f}%)")
    print(f"Variable columns: {variable} ({100*variable/n_col:.2f}%)")
    print(f"Parsimony-informative columns: {pinformative}")
    print(f"Distinct site patterns: {len(site_patterns)}")
    print(f"Gap/missing characters: {gap_pct_overall:.1f}% overall (per-sequence range {gap_range[0]:.1f}%-{gap_range[1]:.1f}%)")
    print(f"Ambiguous IUPAC codes: {ambig_pct_overall:.2f}% of aligned positions\n")

    # --- pairwise identity + haplotype dedup ---
    import random
    pairs = list(itertools.combinations(labels, 2))
    if args.max_pairs and len(pairs) > args.max_pairs:
        random.seed(args.seed)
        pairs = random.sample(pairs, args.max_pairs)
        print(f"(sampling {args.max_pairs} of {len(list(itertools.combinations(labels,2)))} possible pairs for speed)\n")

    idents = []
    for a, b in pairs:
        v = pct_identity(seqs[a], seqs[b])
        if v is not None:
            idents.append(v)
    idents.sort()
    n = len(idents)
    print(f"Pairwise identity across {n} pairs: mean={sum(idents)/n:.1f}%  "
          f"median={idents[n//2]:.1f}%  min={idents[0]:.1f}%  max={idents[-1]:.1f}%")

    unique_seqs = set(seqs.values())
    print(f"Unique haplotypes: {len(unique_seqs)} of {n_seq} "
          f"({100*len(unique_seqs)/n_seq:.1f}%) -- {n_seq - len(unique_seqs)} exact duplicates")


if __name__ == "__main__":
    main()
