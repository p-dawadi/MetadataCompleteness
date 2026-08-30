#!/usr/bin/env python3
"""
parsimony_permutation_test.py

Topology-based Fitch parsimony/permutation test for association between
country of origin and phylogenetic tree topology. Replaces the
adjacency-based tip-order proxy used in the original manuscript
submission (invalid because tip adjacency in a rendered tree depends on
node rotation/ladderization/rendering settings, not only on topology).

Method
------
1. Load the ML tree (../phylogenetics/final_tree/Aligned.Strepto.aln.treefile;
   IQ-TREE 2.4.0, K2P+R2).
2. Match tips to the curated metadata (../data/curated_metadata.csv) and
   assign each a country (see country_mapping.py). Tips without a known
   country are excluded (not treated as an ambiguous/wildcard state).
3. Prune the tree to just the known-country tips (topology only --
   branch lengths are not used in the statistic).
4. Compute the Fitch parsimony score (minimum number of state changes)
   for "country" (a k-state categorical character) on the observed tree,
   using a sequential-pairwise implementation of Fitch's algorithm
   (handles polytomies).
5. Permutation test: shuffle the country labels across the same tips
   (topology held fixed) for N replicates; recompute the parsimony score
   each time. Empirical p = P(null score <= observed score). A lower
   observed score than the null distribution indicates that tips sharing
   a country are more closely related on the tree than chance predicts.

Usage
-----
    python parsimony_permutation_test.py                     # full dataset
    python parsimony_permutation_test.py --exclude-composition-failing
    python parsimony_permutation_test.py --leave-major-study-out
    python parsimony_permutation_test.py --directly-reported-only
    python parsimony_permutation_test.py --n-perm 9999 --seed 20260827

Requires: pandas, dendropy (pip install pandas dendropy)
"""
import argparse
import random
import re
import sys

import pandas as pd
import dendropy

from country_mapping import load_curated_with_country

sys.setrecursionlimit(10000)

TREEFILE = "../phylogenetics/final_tree/Aligned.Strepto.aln.treefile"

# The four sequences that failed IQ-TREE's composition-homogeneity test
# (all whole-genome-shotgun contig-derived records; see manuscript
# Methods: Phylogenetic reconstruction, and docs/data_dictionary.md).
COMPOSITION_FAILING = {
    "NZ_JAQMKN010000038", "QSDO01000081", "VYVX01000034", "PVSY01000043",
}


def base_accession(tip_label):
    """Strip the GenBank version suffix (e.g. 'AF104679.1' -> 'AF104679')
    so tree tip labels can be matched to the metadata table."""
    return re.sub(r"\.\d+$", "", tip_label)


def fitch_score(tree):
    """Sequential-pairwise Fitch algorithm. Requires `.state` to already
    be set on every leaf's taxon. Returns the total parsimony score
    (number of state changes)."""
    cost = [0]

    def downpass(node):
        if node.is_leaf():
            return {node.taxon.state}
        child_sets = [downpass(c) for c in node.child_nodes()]
        current = set(child_sets[0])
        for s in child_sets[1:]:
            intersection = current & s
            if intersection:
                current = intersection
            else:
                current = current | s
                cost[0] += 1
        return current

    downpass(tree.seed_node)
    return cost[0]


def get_largest_study_accessions(df, top_n_countries=4):
    """For each of the `top_n_countries` most-represented countries,
    return the set of accessions belonging to that country's single
    largest contributing "Article name" (used for the
    leave-major-study-out sensitivity analysis)."""
    counts = df["country"].value_counts()
    top_countries = counts.head(top_n_countries).index
    exclude = set()
    for country in top_countries:
        sub = df[df["country"] == country]
        top_article = sub["Article name"].value_counts().index[0]
        exclude |= set(sub[sub["Article name"] == top_article]["Accession Number"])
    return exclude


def run_test(known_labels, acc_to_country, n_perm, seed, label=""):
    tree = dendropy.Tree.get(path=TREEFILE, schema="newick", preserve_underscores=True)
    pruned = tree.clone(depth=1)
    pruned.retain_taxa_with_labels(known_labels)
    pruned.suppress_unifurcations()

    true_states = {}
    for leaf in pruned.leaf_node_iter():
        country = acc_to_country[base_accession(leaf.taxon.label)]
        true_states[leaf.taxon.label] = country
        leaf.taxon.state = country

    observed = fitch_score(pruned)
    n_tips = len(true_states)
    n_states = len(set(true_states.values()))

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
    p_report = f"< {1 / (n_perm + 1):.5f}" if p_value == 0 else f"{p_value:.5f}"

    print(f"--- {label or 'Result'} ---")
    print(f"n_tips={n_tips}  n_states={n_states}  observed_score={observed}")
    print(f"null: mean={null_scores.mean():.2f} sd={null_scores.std():.2f} "
          f"min={null_scores.min()} max={null_scores.max()}  (n_perm={n_perm})")
    print(f"empirical p-value: {p_report}")
    print(f"observed / null_mean = {100 * observed / null_scores.mean():.1f}%\n")

    return {
        "n_tips": n_tips, "n_states": n_states, "observed_score": observed,
        "null_mean": null_scores.mean(), "null_sd": null_scores.std(),
        "null_min": int(null_scores.min()), "null_max": int(null_scores.max()),
        "n_perm": n_perm, "p_value": p_report,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-perm", type=int, default=9999)
    ap.add_argument("--seed", type=int, default=20260827)
    ap.add_argument("--exclude-composition-failing", action="store_true",
                     help="Exclude the 4 sequences that failed the composition-homogeneity test.")
    ap.add_argument("--leave-major-study-out", action="store_true",
                     help="Exclude the single largest contributing study for each of the "
                          "4 most-represented countries (128 records).")
    ap.add_argument("--directly-reported-only", action="store_true",
                     help="Exclude institution-inferred geographic labels; keep only "
                          "directly reported countries.")
    args = ap.parse_args()

    df = load_curated_with_country()
    acc_to_country = dict(zip(df["Accession Number"].astype(str).str.strip(), df["country"]))

    tree = dendropy.Tree.get(path=TREEFILE, schema="newick", preserve_underscores=True)
    all_labels = [l.taxon.label for l in tree.leaf_node_iter()]

    exclude_accessions = set()
    label = "Full dataset"
    if args.exclude_composition_failing:
        exclude_accessions |= COMPOSITION_FAILING
        label = "Excluding composition-test-failing sequences"
    if args.leave_major_study_out:
        exclude_accessions |= get_largest_study_accessions(df)
        label = "Leave-major-study-out"
    if args.directly_reported_only:
        inst_accessions = set(df[df["is_institution_inferred"]]["Accession Number"])
        exclude_accessions |= inst_accessions
        label = "Directly-reported locations only"

    known_labels = [
        lbl for lbl in all_labels
        if acc_to_country.get(base_accession(lbl)) is not None
        and base_accession(lbl) not in exclude_accessions
    ]

    run_test(known_labels, acc_to_country, args.n_perm, args.seed, label=label)


if __name__ == "__main__":
    main()
