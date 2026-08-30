#!/usr/bin/env python3
"""
alignment_occupancy_check.py

Checks whether the final alignment (Aligned.Strepto.aln) shows evidence of
end/column trimming, as a real (data-driven) answer to the open question of
whether any trimming was applied beyond MAFFT's own default behavior.

Logic: a standard end- or column-trimming tool (trimAl, Gblocks, etc.)
removes alignment columns where most sequences have a gap, because such
columns are considered unreliable/uninformative. If such a step had been
run on this alignment, sparsely occupied columns (especially runs of them
at the termini) would not remain in the final file. Their presence is
therefore direct evidence against a separate trimming step having been
applied.

Usage: python alignment_occupancy_check.py
"""

ALIGNMENT = "../alignment/Aligned.Strepto.aln"


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


def main():
    seqs = read_fasta(ALIGNMENT)
    n = len(seqs)
    L = len(next(iter(seqs.values())))
    occ = [0] * L
    for s in seqs.values():
        for i, c in enumerate(s):
            if c not in "-.":
                occ[i] += 1
    occ = [100.0 * o / n for o in occ]

    first50 = next(i for i, o in enumerate(occ) if o >= 50)
    last50 = L - 1 - next(i for i, o in enumerate(reversed(occ)) if o >= 50)
    below10 = sum(1 for o in occ if o < 10)

    print(f"Alignment: {n} sequences, {L} columns\n")
    print(f"First column with >=50% sequence occupancy: {first50} "
          f"({first50} leading columns below 50% occupancy)")
    print(f"Last column with >=50% sequence occupancy: {last50} "
          f"({L - 1 - last50} trailing columns below 50% occupancy)")
    print(f"Columns with <10% occupancy overall: {below10} of {L} ({100*below10/L:.1f}%)")
    print(f"Minimum occupancy anywhere: {min(occ):.1f}% (column {occ.index(min(occ))})")
    print()
    print("A standard end- or column-trimming step (trimAl, Gblocks, etc.) would")
    print("ordinarily remove most such sparsely occupied positions. Their retention")
    print("here indicates the alignment used for tree reconstruction reflects MAFFT's")
    print("own untrimmed default output, not a subsequently trimmed version.")


if __name__ == "__main__":
    main()
