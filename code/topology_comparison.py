#!/usr/bin/env python3
"""
topology_comparison.py

Quantifies how much tree topology actually changes between the primary
tree (all 350 curated records) and the sensitivity tree (346 records,
excluding the 4 identity/composition-outlier WGS-contig sequences),
using the Robinson-Foulds (RF) distance -- the standard measure of
topological difference between two trees (count of bipartitions/splits
present in one tree but not the other).

Since the two trees don't have identical taxon sets (350 vs 346 tips),
the primary tree is first pruned down to the 346 shared tips before
computing RF distance, so the comparison is like-for-like (topology only,
not affected by the 4 tips' mere presence/absence).

Reports:
  - raw RF distance (number of differing bipartitions)
  - normalized RF distance (0 = identical topology, 1 = maximally different)
  - the maximum possible RF distance for a tree this size (for context)
  - OPTIONALLY (--min-support), the same comparison after first collapsing
    every internal branch below a UFBoot support threshold on each tree
    (independently, using each tree's own support values, before pruning
    to shared taxa). Raw RF distance treats a disagreement among
    essentially unresolved/near-polytomy branches exactly the same as a
    disagreement among rock-solid ones; collapsing low-support branches
    first answers the more relevant question of whether the *confidently
    resolved* backbone topology agrees between the two trees.

Usage:
    python topology_comparison.py --primary ../../tree_350_primary.treefile \\
        --sensitivity ../../tree_346_sensitivity.treefile
    python topology_comparison.py --primary ../../tree_350_primary.treefile \\
        --sensitivity ../../tree_346_sensitivity.treefile --min-support 70
"""
import argparse
import dendropy
from dendropy.calculate import treecompare


def collapse_low_support(tree, min_support):
    """Collapse every internal branch (edge) whose UFBoot support label is
    below min_support, merging its children up into its parent (creating a
    polytomy there). Must be called on the tree BEFORE any taxon pruning --
    support values are attached to the original tree's internal nodes, and
    pruning/suppressing unifurcations can merge nodes in ways that make a
    post-hoc support lookup ambiguous. Leaves and the root are never
    touched. Internal nodes with no support label (e.g. unlabeled
    zero-length branches within a cluster of identical sequences) are left
    as-is rather than being treated as automatically low-confidence."""
    to_collapse = []
    for node in tree.preorder_node_iter():
        if node.parent_node is None or node.is_leaf():
            continue
        label = node.label
        if label is None or str(label).strip() == "":
            continue
        try:
            support = float(label)
        except ValueError:
            continue
        if support < min_support:
            to_collapse.append(node)
    for node in to_collapse:
        parent = node.parent_node
        if parent is None:
            continue
        for child in list(node.child_nodes()):
            node.remove_child(child)
            parent.add_child(child)
        parent.remove_child(node)
    tree.encode_bipartitions()


def compare(t_primary, t_sens, shared, label):
    t_primary_pruned = t_primary.clone(depth=1)
    t_primary_pruned.retain_taxa_with_labels(shared)
    t_primary_pruned.suppress_unifurcations()
    t_sens_pruned = t_sens.clone(depth=1)
    t_sens_pruned.retain_taxa_with_labels(shared)
    t_sens_pruned.suppress_unifurcations()

    t_primary_pruned.encode_bipartitions()
    t_sens_pruned.encode_bipartitions()

    rf = treecompare.symmetric_difference(t_primary_pruned, t_sens_pruned)
    n_tips = len(shared)
    max_rf = 2 * (n_tips - 3)
    norm_rf = rf / max_rf if max_rf > 0 else 0.0

    print(f"--- {label} ---")
    print(f"Robinson-Foulds distance (on {n_tips} shared tips): {rf}")
    print(f"Maximum possible RF distance for {n_tips} tips: {max_rf}")
    print(f"Normalized RF distance: {norm_rf:.4f} "
          f"(0 = identical topology, 1 = maximally different)\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--primary", required=True, help="Path to tree_350_primary.treefile")
    ap.add_argument("--sensitivity", required=True, help="Path to tree_346_sensitivity.treefile")
    ap.add_argument("--min-support", type=float, default=None,
                     help="If set, also report RF distance after collapsing internal branches "
                          "with UFBoot support below this threshold (e.g. 70) on each tree first.")
    args = ap.parse_args()

    tns = dendropy.TaxonNamespace()
    t_primary = dendropy.Tree.get(path=args.primary, schema="newick",
                                   taxon_namespace=tns, preserve_underscores=True)
    t_sens = dendropy.Tree.get(path=args.sensitivity, schema="newick",
                                taxon_namespace=tns, preserve_underscores=True)

    sens_labels = {leaf.taxon.label for leaf in t_sens.leaf_node_iter()}
    primary_labels = {leaf.taxon.label for leaf in t_primary.leaf_node_iter()}
    shared = sens_labels & primary_labels
    only_primary = primary_labels - sens_labels

    print(f"Primary tree: {len(primary_labels)} tips")
    print(f"Sensitivity tree: {len(sens_labels)} tips")
    print(f"Shared tips: {len(shared)}")
    print(f"Tips only in primary (expected: the 4 excluded outliers): {sorted(only_primary)}\n")

    compare(t_primary, t_sens, shared, "Full topology (all resolved branches)")

    if args.min_support is not None:
        t_primary2 = t_primary.clone(depth=1)
        t_sens2 = t_sens.clone(depth=1)
        collapse_low_support(t_primary2, args.min_support)
        collapse_low_support(t_sens2, args.min_support)
        compare(t_primary2, t_sens2, shared,
                f"Support-filtered (branches with UFBoot < {args.min_support}% collapsed first)")


if __name__ == "__main__":
    main()
