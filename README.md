# *Streptococcus anginosus* 16S rRNA Metadata Completeness & Phylogenetic Structure — Analysis Pipeline

This repository contains the analysis pipeline, curated data, alignments, tree files, and
supporting scripts behind the manuscript *"Metadata Completeness and Phylogenetic Structure
of Publicly Available Streptococcus anginosus 16S rRNA Gene Sequences"* (BMC Microbiology,
under review). It does **not** include the manuscript text itself — only the underlying
data and code needed to reproduce and verify the reported analyses.

## Repository structure

```
code/                            Analysis scripts (run from inside this folder)
  1_retrieve_sequences.py          Step 1: NCBI retrieval
  4_alignment_stats.py             Step 4: alignment statistics
  5_alignment_occupancy_check.py   Step 5: alignment trimming/occupancy check
  country_mapping.py               Geographic-location -> country standardization (imported by other scripts)
  composition_anomaly_check.py     Network-free composition/assembly QC on the 4 flagged records
  check_no_reference_sequences.py  Confirms no external/reference sequences leaked into the tip set
  parsimony_permutation_test.py    Fitch parsimony / permutation test (country vs. topology)
  phylo_uncertainty_check.py       Same test repeated across all UFBoot trees
  topology_comparison.py           Robinson-Foulds distance between the two trees
  bootstrap_support_distribution.py  UFBoot support-value summary

data/
  curated_metadata.csv             Curated metadata table (350 records)
  sequences/
    Strepto_350_curated_only.fasta   All 350 curated records
    Strepto_346_excl_outliers.fasta  346 records (350 minus the 4 composition-failing WGS contigs)
  alignments/
    Aligned_350_primary.fasta        MAFFT alignment, 350-record dataset
    Aligned_346_sensitivity.fasta    MAFFT alignment, 346-record dataset
  silva_screening/
    silva_hits.tsv                   BLAST hits vs. the SILVA SSU rRNA reference database
    silva_best_hits.tsv              Best hit per query

trees/
  350_primary/                     IQ-TREE output, primary (350-record) analysis
  346_sensitivity/                 IQ-TREE output, sensitivity (346-record) analysis, plus
                                    phylo_uncertainty_346.log (bootstrap-tree uncertainty result)
```

## Workflow

### 1. Sequence retrieval
`code/1_retrieve_sequences.py` retrieves *Streptococcus anginosus* 16S rRNA records from
NCBI Nucleotide via Biopython's Entrez module (one record at a time, with retry/resume
logic for unreliable network conditions), restricted to sequences 1,200–1,700 bp in
length, then tabulates results with pandas.

```
python code/1_retrieve_sequences.py
```

Requires `biopython`, `pandas`, `requests`. After deduplication and removal of blank
export artifacts, 350 unique records remained for analysis (lengths 1,207–1,573 bp).

> **Note on `retrieval.log`:** the original retrieval run's console log was not preserved
> on disk; only a later, failed re-verification attempt (a missing-dependency error) was
> found when this repository was assembled, so it has been left out rather than included
> as evidence of the original run. The retrieval script, query logic, and length filters
> above are unchanged from what actually produced the 350-record dataset.

### 2. Metadata curation
Host, geographic origin, and sample source were extracted from each GenBank record (and,
where available, its associated publication) and standardized into `data/curated_metadata.csv`.
Geographic entries naming an institution rather than a country are mapped to a country by
`code/country_mapping.py`, which flags such rows via `is_institution_inferred` so they can
be excluded in the "directly-reported-only" sensitivity analysis (see step 6).

### 3. Taxonomic/composition screening
Sequences were screened against the SILVA SSU rRNA reference database (BLAST) to check
for mislabeled or non-target sequences; results are in `data/silva_screening/`.
`code/composition_anomaly_check.py` runs additional network-free composition/assembly
checks (GC content, ambiguous-base content, repeat content) on the four records that
failed IQ-TREE's composition-homogeneity test:

```
python code/composition_anomaly_check.py data/sequences/Strepto_350_curated_only.fasta
```

`code/check_no_reference_sequences.py` confirms every tip label in an alignment/tree
traces back to a real curated accession, with no external reference sequence left in by
accident:

```
python code/check_no_reference_sequences.py data/alignments/Aligned_350_primary.fasta
```

### 4. Sequence alignment
The curated sequences were aligned with MAFFT (command-line) to produce
`data/alignments/Aligned_350_primary.fasta` (350 records) and
`data/alignments/Aligned_346_sensitivity.fasta` (346 records, excluding the four
composition-failing whole-genome-shotgun contig-derived records
NZ_JAQMKN010000038, QSDO01000081, VYVX01000034, PVSY01000043).

### 5. Alignment QC
`code/4_alignment_stats.py` computes alignment statistics (column count, constant/
parsimony-informative sites, gap and ambiguous-base content, pairwise identity, unique
haplotypes) directly from the MAFFT output:

```
python code/4_alignment_stats.py data/alignments/Aligned_350_primary.fasta
python code/4_alignment_stats.py data/alignments/Aligned_346_sensitivity.fasta
```

`code/5_alignment_occupancy_check.py` checks whether the alignment shows evidence of
end/column trimming beyond MAFFT's own default output:

```
python code/5_alignment_occupancy_check.py data/alignments/Aligned_350_primary.fasta
```

### 6. Phylogenetic reconstruction
Maximum-likelihood trees were built with **IQ-TREE 3.1.3**, one per dataset. The exact
commands, taken directly from each run's own recorded log, were:

```
# 350-record primary tree — ModelFinder selects the best-fit model (BIC: K3Pu+F+R3)
iqtree -s data/alignments/Aligned_350_primary.fasta -st DNA -m MFP -bb 1000 -nt AUTO -seed 20260827 -pre tree_350_primary

# 346-record sensitivity tree — TPM3u+F+R3, with UFBoot trees written out for step 8
iqtree -s data/alignments/Aligned_346_sensitivity.fasta -st DNA -m TPM3u+F+R3 -bb 1000 -wbt -nt AUTO -seed 20260827 -pre tree_346_sensitivity
```

Both runs used 1,000 ultrafast bootstrap replicates. Full output (`.iqtree` report,
`.log`, `.treefile`, `.contree`, `.bionj`, `.mldist`, `.model.gz`, `.splits.nex`,
`.uniqueseq.phy`) is in `trees/350_primary/` and `trees/346_sensitivity/`;
`trees/346_sensitivity/` additionally has `tree_346_sensitivity.ufboot` (the 1,000
bootstrap trees, needed for step 8). IQ-TREE checkpoint files (`.ckp.gz`) are excluded —
they only support resuming an interrupted run and carry no analytical content.

> **Note on `run_iqtree.pbs`:** an older PBS job script from an earlier iteration of this
> analysis (different model/settings, different input filename) is not included here, to
> avoid misrepresenting how the final trees above were actually produced. The commands
> above are copied verbatim from the "Command:" line each `.log` file itself records.

### 7. Country / topology association (Fitch parsimony + permutation test)
`code/parsimony_permutation_test.py` tests whether tips sharing a reported country are
more closely related on the tree than expected by chance: it computes the Fitch
parsimony score for "country" as a categorical character on the tree, then compares it
to a null distribution built by permuting country labels across tips (topology fixed).

```
python code/parsimony_permutation_test.py --treefile trees/350_primary/tree_350_primary.treefile
python code/parsimony_permutation_test.py --treefile trees/346_sensitivity/tree_346_sensitivity.treefile

# sensitivity variants
python code/parsimony_permutation_test.py --treefile trees/350_primary/tree_350_primary.treefile --exclude-composition-failing
python code/parsimony_permutation_test.py --treefile trees/350_primary/tree_350_primary.treefile --leave-major-study-out
python code/parsimony_permutation_test.py --treefile trees/350_primary/tree_350_primary.treefile --directly-reported-only
```

Requires `pandas`, `dendropy`. Reads tip-to-country assignments from
`data/curated_metadata.csv` via `country_mapping.py`.

### 8. Phylogenetic uncertainty across bootstrap trees
`code/phylo_uncertainty_check.py` repeats the identical Fitch parsimony/permutation test
across all 1,000 UFBoot trees (not just the single best ML tree), so the
observed/null-ratio and p-value are reported as a distribution rather than a single-tree
point estimate — directly addressing how sensitive the topology-association result is to
the point-estimate tree:

```
python code/phylo_uncertainty_check.py --ufboot trees/346_sensitivity/tree_346_sensitivity.ufboot
```

The result of this run is recorded in `trees/346_sensitivity/phylo_uncertainty_346.log`.

### 9. Topology comparison (Robinson-Foulds distance)
`code/topology_comparison.py` quantifies how much the primary (350-tip) and sensitivity
(346-tip) trees actually differ in topology, using the Robinson-Foulds distance after
pruning the primary tree to the 346 shared tips:

```
python code/topology_comparison.py --primary trees/350_primary/tree_350_primary.treefile --sensitivity trees/346_sensitivity/tree_346_sensitivity.treefile
python code/topology_comparison.py --primary trees/350_primary/tree_350_primary.treefile --sensitivity trees/346_sensitivity/tree_346_sensitivity.treefile --min-support 70
```

Requires `dendropy`.

### 10. Bootstrap support distribution
`code/bootstrap_support_distribution.py` summarizes the UFBoot support values on the
internal branches of a treefile:

```
python code/bootstrap_support_distribution.py trees/350_primary/tree_350_primary.treefile
python code/bootstrap_support_distribution.py trees/346_sensitivity/tree_346_sensitivity.treefile
```

## Dependencies

- Python 3 with `biopython`, `pandas`, `dendropy`
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) (command-line)
- [IQ-TREE](http://www.iqtree.org/) 3.1.3
- BLAST+ (for the SILVA screening step)

## Citation

If you use this data or code, please cite:

> Pokharel B, Nepal S, Roy M, Regmi S, Dawadi P. (2026).
> **Metadata Completeness and Phylogenetic Structure of Publicly Available *Streptococcus
> anginosus* 16S rRNA Gene Sequences.**
> *Short title: Metadata gaps and phylogenetic structure in public S. anginosus 16S archives.*
> Manuscript under review.

A full citation with journal, volume, and DOI will be added here once the manuscript is
published.

## Maintainer

Pipeline developed and maintained by Prabin Dawadi (2026).

This repository accompanies a manuscript currently under peer review; a formal open-data
license (e.g., CC BY 4.0 for data, MIT for code) will be applied here once the journal's
publication and licensing terms are finalized.
