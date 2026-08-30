#!/usr/bin/env python3
"""
composition_anomaly_check.py

A further, real (network-free) check on the four sequences that failed
IQ-TREE's composition-homogeneity test (NZ_JAQMKN010000038.1,
QSDO01000081.1, VYVX01000034.1, PVSY01000043.1). Formal BLAST/CheckM
screening against external databases was not available in this
environment; this script instead computes, directly from the sequences
already in hand, several composition/assembly-quality signals that would
flag misassembly or contamination if present:

  - GC content (compared with the dataset-wide distribution)
  - ambiguous/N-base content
  - ungapped sequence length
  - the longest exact repeated substring >= 20 bp within each sequence
    (a simple, conservative signal for assembly duplication/misjoin --
    real biological repeats exist too, so this flags candidates for
    review rather than proving misassembly)

Usage: python composition_anomaly_check.py
"""
import re
from collections import Counter

ALIGNMENT = "../alignment/Aligned.Strepto.aln"
FLAGGED = {"NZ_JAQMKN010000038.1", "QSDO01000081.1", "VYVX01000034.1", "PVSY01000043.1"}


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


def ungapped(seq):
    return seq.replace("-", "").replace(".", "")


def gc_content(seq):
    s = seq.upper()
    gc = sum(1 for c in s if c in "GC")
    acgt = sum(1 for c in s if c in "ACGT")
    return 100.0 * gc / acgt if acgt else float("nan")


def n_content(seq):
    s = seq.upper()
    n = sum(1 for c in s if c not in "ACGT")
    return 100.0 * n / len(s) if s else float("nan")


def longest_exact_repeat(seq, min_len=20, max_check=40):
    """Longest exact repeated substring of length >= min_len, found by
    a simple expanding search (adequate for sequences of this size;
    not a suffix-array implementation)."""
    s = seq.upper()
    best = 0
    seen = {}
    k = min_len
    while k <= max_check:
        kmers = Counter(s[i:i + k] for i in range(len(s) - k + 1))
        if any(v > 1 for v in kmers.values()):
            best = k
            k += 5
        else:
            break
    return best


def main():
    seqs = read_fasta(ALIGNMENT)
    stats = {}
    for label, seq in seqs.items():
        u = ungapped(seq)
        stats[label] = {
            "len": len(u),
            "gc": gc_content(u),
            "n_pct": n_content(u),
            "longest_repeat": longest_exact_repeat(u),
        }

    all_gc = [v["gc"] for v in stats.values()]
    all_n = [v["n_pct"] for v in stats.values()]
    mean_gc, sd_gc = sum(all_gc) / len(all_gc), (sum((x - sum(all_gc) / len(all_gc))**2 for x in all_gc) / len(all_gc)) ** 0.5
    mean_n = sum(all_n) / len(all_n)

    print(f"Dataset-wide (n={len(stats)}): mean GC = {mean_gc:.2f}% (SD {sd_gc:.2f}); mean ambiguous-base content = {mean_n:.3f}%\n")
    print("Flagged (composition-test-failing) sequences:")
    print(f"{'accession':22s} {'len':>6s} {'GC%':>7s} {'GC z-score':>11s} {'ambig%':>8s} {'longest exact repeat (bp)':>26s}")
    for label in sorted(FLAGGED):
        if label not in stats:
            print(f"  [not found in alignment: {label}]")
            continue
        s = stats[label]
        z = (s["gc"] - mean_gc) / sd_gc if sd_gc else float("nan")
        print(f"{label:22s} {s['len']:6d} {s['gc']:7.2f} {z:11.2f} {s['n_pct']:8.3f} {s['longest_repeat']:26d}")

    print("\nFor comparison, dataset-wide GC z-score distribution (all 354 sequences):")
    zs = [(v["gc"] - mean_gc) / sd_gc for v in stats.values() if sd_gc]
    zs_sorted = sorted(zs)
    print(f"  min z = {zs_sorted[0]:.2f}, max z = {zs_sorted[-1]:.2f}, "
          f"n with |z| > 2 = {sum(1 for z in zs if abs(z) > 2)} of {len(zs)}")


if __name__ == "__main__":
    main()
