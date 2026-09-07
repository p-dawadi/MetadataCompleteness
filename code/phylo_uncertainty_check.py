#!/usr/bin/env python3
"""
phylo_uncertainty_check.py

Quantifies how much the country/topology association result (see
parsimony_permutation_test.py) depends on the single point-estimate ML
tree, by re-running the identical Fitch parsimony/permutation test across
the full set of ultrafast bootstrap (UFBoot) trees from the same IQ-TREE
run (the `.ufboot` file IQ-TREE writes automatically when run with -bb).

This directly answers the editor's comment that "phylogenetic uncertainty
... remain[s] incompletely addressed": rather than reporting a single
observed/null ratio and p-value from one tree, this reports their
distribution across many plausible trees consistent with the data's
bootstrap support, so a reader can see whether the association is a
robust signal or an artifact of one particular tree estimate.

Method
------
Identical to parsimony_permutation_test.py's main test (same Fitch
algorithm, same country assignment, same permutation logic), but looped
over every tree in the .ufboot file instead of the single best ML tree.
A smaller per-tree permutation count is used (default 999) to keep total
runtime reasonable across ~1000 trees; this is still enough to place the
observed score in the null distribution at coarse (0.1%) resolution per
tree, and the *aggregate* result (the distribution of ratios/p-values
across all trees) is the quantity of interest here, not any single tree's
p-value precision.

Usage
-----
    python phylo_uncertainty_check.py --ufboot ../../tree_350_primary.ufboot
    python phylo_uncertainty_check.py --ufboot ../../tree_350_primary.ufboot \\
        --exclude-composition-failing --n-trees 1000 --n-perm 999

Requires: pandas, dendropy (same as parsimony_permutation_test.py)
"""
import argparse
import random
import re
import sys

import pandas as pd
import dendropy

from country_mapping import load_curated_with_country
from parsimony_permutation_test import fitch_score, base_accession, COMPOSITION_FAILING

sys.setrecursionlimit(10000)


def run_one_tree(tree, known_labels, acc_to_country, n_perm, seed):
    pruned = tree.clone(depth=1)
    pruned.retain_taxa_with_labels(known_labels)
    pruned.suppress_unifurcations()

    true_states = {}
    for leaf in pruned.leaf_node_iter():
        country = acc_to_country[base_accession(leaf.taxon.label)]
        true_states[leaf.taxon.label] = country
        leaf.taxon.state = country

    observed = fitch_score(pruned)

    random.seed(seed)
    labels_list = list(true_states.keys())
    states_list = [true_states[l] for l in labels_list]
    leaf_nodes = list(pruned.leaf_node_iter())
    null_scores = []
    for _ in range(n_perm):
        shuffled = states_list[:]
        random.shuffle(shuffled)
        perm_map = dict(zip(labels_list, shuffled))
        for leaf in leaf_nodes:
            leaf.taxon.state = perm_map[leaf.taxon.label]
        null_scores.append(fitch_score(pruned))
    null_scores = pd.Series(null_scores)

    p_value = (null_scores <= observed).sum() / len(null_scores)
    ratio = 100 * observed / null_scores.mean()
    return observed, null_scores.mean(), p_value, ratio


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ufboot", required=True, help="Path to IQ-TREE's .ufboot multi-tree file")
    ap.add_argument("--n-trees", type=int, default=None,
                     help="Use only the first N trees from the .ufboot file (default: all)")
    ap.add_argument("--n-perm", type=int, default=999,
                     help="Permutations per tree (default 999; kept lower than the main "
                          "test's 9999 since this runs across ~1000 trees)")
    ap.add_argument("--seed", type=int, default=20260827)
    ap.add_argument("--exclude-composition-failing", action="store_true")
    ap.add_argument("--sig-threshold", type=float, default=0.05)
    args = ap.parse_args()

    df = load_curated_with_country()
    acc_to_country = dict(zip(df["Accession Number"].astype(str).str.strip(), df["country"]))

    print(f"Loading bootstrap trees from {args.ufboot} ...")
    trees = dendropy.TreeList.get(path=args.ufboot, schema="newick", preserve_underscores=True)
    if args.n_trees:
        trees = trees[:args.n_trees]
    print(f"Loaded {len(trees)} trees.\n")

    exclude_accessions = COMPOSITION_FAILING if args.exclude_composition_failing else set()

    all_labels = [l.taxon.label for l in trees[0].leaf_node_iter()]
    known_labels = [
        lbl for lbl in all_labels
        if acc_to_country.get(base_accession(lbl)) is not None
        and base_accession(lbl) not in exclude_accessions
    ]
    print(f"n_tips used per tree: {len(known_labels)}\n")

    ratios, pvals, observeds = [], [], []
    for i, tree in enumerate(trees):
        observed, null_mean, p_value, ratio = run_one_tree(
            tree, known_labels, acc_to_country, args.n_perm, args.seed + i
        )
        observeds.append(observed)
        ratios.append(ratio)
        pvals.append(p_value)
        if (i + 1) % 100 == 0 or i == len(trees) - 1:
            print(f"  ...{i+1}/{len(trees)} trees done", file=sys.stderr)

    ratios = pd.Series(ratios)
    pvals = pd.Series(pvals)
    n_sig = (pvals < args.sig_threshold).sum()

    print("=== Phylogenetic-uncertainty summary across bootstrap trees ===")
    print(f"Trees tested: {len(trees)}  (permutations per tree: {args.n_perm})")
    print(f"observed/null_mean %% ratio: median={ratios.median():.1f}%%  "
          f"IQR=[{ratios.quantile(0.25):.1f}%%, {ratios.quantile(0.75):.1f}%%]  "
          f"range=[{ratios.min():.1f}%%, {ratios.max():.1f}%%]")
    print(f"empirical p-value: median={pvals.median():.4f}  "
          f"max={pvals.max():.4f}")
    print(f"trees with p < {args.sig_threshold}: {n_sig}/{len(trees)} ({100*n_sig/len(trees):.1f}%%)")


if __name__ == "__main__":
    main()
