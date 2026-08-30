#!/usr/bin/env python3
"""
Summarizes the ultrafast-bootstrap (UFBoot) support distribution across the
internal branches of the final maximum-likelihood tree
(phylogenetics/final_tree/Aligned.Strepto.aln.treefile), computed directly
from the real IQ-TREE 2.4.0 run (iqtree2 -s Aligned.Strepto.aln -m K2P+R2
-nt 8 -bb 10000 -redo -safe; see phylogenetics/final_tree/run_iqtree.pbs).

Support values are written by IQ-TREE as an integer immediately after the
closing parenthesis of an internal node, before its branch length
(e.g. ")99:0.0015128513"). Not every internal branch carries a printed
support value: branches with essentially zero length within a cluster of
identical or near-identical sequences represent an arbitrary resolution of
what is effectively a polytomy, and IQ-TREE does not assign a separate
support value to these. Both counts are reported below rather than treating
the unlabeled branches as 0% or excluding them silently.

Run: python bootstrap_support_distribution.py
"""
import re
import statistics
import pathlib

TREEFILE = pathlib.Path(__file__).resolve().parent.parent / "phylogenetics" / "final_tree" / "Aligned.Strepto.aln.treefile"

def main():
    tree = TREEFILE.read_text().strip()

    after_paren = re.findall(r"\)([^,()]*)", tree)
    total_internal_nodes = len(after_paren)  # includes the basal trifurcation point of the unrooted tree

    supports = []
    unlabeled = 0
    for tok in after_paren:
        m = re.match(r"^(\d+):", tok)
        if m:
            supports.append(int(m.group(1)))
        elif re.match(r"^:", tok):
            unlabeled += 1
        # else: the final ';' terminator of the newick string (the root point itself) -- not a real branch

    n = len(supports)
    print(f"Total internal nodes parsed from treefile: {total_internal_nodes}")
    print(f"  Real internal branches (excluding the basal trifurcation point): {n + unlabeled}")
    print(f"  With an explicit UFBoot support value: {n} ({100*n/(n+unlabeled):.1f}%)")
    print(f"  Zero-length, unlabeled (polytomy-like) branches: {unlabeled} ({100*unlabeled/(n+unlabeled):.1f}%)")
    print()
    print(f"Support value distribution (n={n} labeled branches):")
    print(f"  Range: {min(supports)}-{max(supports)}%")
    print(f"  Mean: {statistics.mean(supports):.1f}%")
    print(f"  Median: {statistics.median(supports)}%")
    print()
    for t in (50, 70, 80, 90, 95, 99, 100):
        k = sum(1 for s in supports if s >= t)
        print(f"  >= {t}% support: {k}/{n} ({100*k/n:.1f}%)")

if __name__ == "__main__":
    main()
