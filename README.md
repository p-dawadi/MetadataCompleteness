# *Streptococcus anginosus* 16S rRNA phylogeography

Data, alignment, phylogenetic tree, analysis code, and manuscript files for a
study of publicly available *Streptococcus anginosus* 16S rRNA gene sequences
(1200–1700 bp) retrieved from NCBI Nucleotide, testing for an association
between geographic origin and phylogenetic position.

[AUTHOR TO CONFIRM: manuscript title, author list, and journal/DOI once
assigned — add a "How to cite" section below once available.]

## Repository structure

```
data/            Curated metadata (350 records), accession list, raw working copy
alignment/       Final multiple sequence alignment (FASTA) used for tree inference
phylogenetics/   IQ-TREE / ModelFinder outputs: final ML tree, run logs, PBS job scripts
code/            Python analysis scripts (country mapping, parsimony/permutation test)
results/         Output of the parsimony/permutation test, all four analysis variants
figures/         Main-text Figures 1-5 (metadata completeness / breakdown charts)
docs/            Data dictionary documenting all fields and derived columns
manuscript/      Manuscript (clean and tracked-changes) and the reviewer response letter
```

### `data/`

- `curated_metadata.xlsx` / `curated_metadata.csv` — the 350 curated records
  (accession, title, sequence length, sample source, NCBI link, geographic
  location, host, source article). See [`docs/data_dictionary.md`](docs/data_dictionary.md)
  for column definitions.
- `accession_list.txt` — the 350 accession numbers, one per line.
- `raw_working_metadata_PRE-PUBLICATION.xlsx` — the original working
  spreadsheet (with curation color-coding), kept for audit-trail purposes.

### `alignment/`

- `Aligned.Strepto.aln` — the final multiple sequence alignment (354
  sequences, 2,310 aligned positions), built de novo with MAFFT v7.526 in
  UGENE desktop (v53.1) among the retrieved sequences themselves. The
  candidate sequences were screened against the SILVA ribosomal RNA gene
  database beforehand to confirm 16S rRNA gene identity (see manuscript
  Methods: Sequence alignment); SILVA was not itself the MAFFT alignment
  reference, and no SILVA reference sequences are present among the 354
  aligned/tree sequences (see `code/check_no_reference_sequences.py`) —
  used as input to IQ-TREE.

### `phylogenetics/`

- `final_tree/` — the maximum-likelihood tree (`Aligned.Strepto.aln.treefile`,
  Newick format, ultrafast-bootstrap support values annotated at internal
  nodes; IQ-TREE 2.4.0, K2P+R2 substitution model, 10,000 ultrafast
  bootstrap replicates — `run_iqtree.pbs` shows the exact command,
  `iqtree -s Aligned.Strepto.aln -m K2P+R2 -nt 8 -bb 10000 -redo -safe`),
  the full IQ-TREE run report (`.iqtree`), the UFBoot consensus tree
  (`.contree`) and split support file (`.splits.nex`), run log, PBS job
  script, and job stdout.
- `modelfinder/` — a separate, earlier exploratory run (`iqtree -m MFP -nt
  AUTO`, no bootstrap) used only for substitution-model comparison across
  AIC/AICc/BIC; not the run used for any reported tree or support value.
  Because this run predates the composition-based sequence review, it
  reports 189 near-zero internal branches rather than the 195 reported in
  the manuscript and in `final_tree/` — both are real IQ-TREE outputs, just
  from two different runs; the manuscript's figure (195) is always the one
  from the final bootstrap run in `final_tree/`, which is the one used
  throughout the analysis.

### `code/`

- `country_mapping.py` — maps the free-text `Geographical location` metadata
  field to a standardized country name, and flags records where the country
  was inferred from a submitting institution rather than stated directly.
  Run directly (`python country_mapping.py`) for a self-test summary.
- `parsimony_permutation_test.py` — the topology-based Fitch parsimony /
  label-permutation test for association between country of origin and tree
  topology (see [`results/parsimony_test_results.md`](results/parsimony_test_results.md)
  for the full method description and results). Run with `--help` for all
  options.
- `type_strain_identity_check.py` — computes ungapped pairwise percent
  identity, directly from the alignment, of every curated sequence against
  the two *S. anginosus* type-strain records already present in the dataset
  (`AF104678.1` / ATCC 33397 and `NR_041722.2` / SK52 = DSM 20563). A
  network-free substitute for BLAST/type-strain verification (see
  [`results/type_strain_identity_check.txt`](results/type_strain_identity_check.txt)).
- `composition_anomaly_check.py` — for the four sequences that failed
  IQ-TREE's composition-homogeneity test, computes GC content (vs. the
  dataset-wide distribution), ambiguous-base content, and the longest exact
  internal repeat, as a further network-free inspection for signs of
  misassembly (see [`results/composition_anomaly_check.txt`](results/composition_anomaly_check.txt)).
- `check_no_reference_sequences.py` — confirms that all 354 tip labels in
  the final alignment/tree match either a curated study accession or one of
  the four known additional GenBank records, with zero unexplained
  sequences — consistent with an alignment built solely from the retrieved
  sequences, with no external reference sequences added (see
  [`results/check_no_reference_sequences.txt`](results/check_no_reference_sequences.txt)).
- `alignment_occupancy_check.py` — computes per-column sequence occupancy
  across the alignment to check for trimming: a trimming tool would remove
  sparsely occupied columns, and hundreds remain (806 of 2,310 columns below
  10% occupancy), indicating the alignment is MAFFT's untrimmed default
  output (see [`results/alignment_occupancy_check.txt`](results/alignment_occupancy_check.txt)).
- `bootstrap_support_distribution.py` — parses the UFBoot support value
  annotated at every internal node of `phylogenetics/final_tree/Aligned.Strepto.aln.treefile`
  and reports the full distribution (range, mean, median) and the
  proportion of branches exceeding common support thresholds (50/70/80/90/95/99/100%),
  distinguishing branches with an explicit support value from zero-length,
  unlabeled branches within near-identical-sequence clusters (see
  [`results/bootstrap_support_distribution.txt`](results/bootstrap_support_distribution.txt)).

### `results/`

- `parsimony_test_results.md` — combined results table and interpretation for
  all four analysis variants (full dataset; excluding the four
  composition-test-failing sequences; leave-major-study-out; directly-reported
  locations only).
- `full_dataset.txt`, `exclude_composition_failing.txt`,
  `leave_major_study_out.txt`, `directly_reported_only.txt` — raw console
  output from each run.

### `figures/`

- `Figure1_metadata_completeness_bar.png` / `Figure2_metadata_completeness_pies.png`
  — reported-vs-missing metadata counts for sample source, geographic
  origin, and host (n=350), as a bar chart and as three pie charts.
- `Figure3_geographic_origin_by_country.png` — geographic origin by country
  among the 274 records with reported location (19 countries).
- `Figure4_host_category.png` — host category among the 292 records with
  reported host metadata (289 human, 1 dog, 1 sea water, 1 wastewater).
- `Figure5_sample_source_by_category.png` — sample source among the 242
  records with reported source metadata, grouped into 9 clinical/anatomical
  categories.

All five were generated by `code/make_metadata_figures.py` directly from
`data/curated_metadata.csv` (and, for Figure 3, `code/country_mapping.py`);
see that script's docstring for exactly which counts were independently
re-derived from the raw data versus taken from the manuscript's own
already-verified totals.

### `docs/`

- `data_dictionary.md` — full column-by-column documentation of the curated
  metadata table, including the country-mapping rules, known limitations,
  and the type-strain species-identification check.
- `institution_to_country_mapping.csv` / `Supplementary_File_2_Data_Dictionary.xlsx`
  — the institution-to-country mapping table and sample-source category
  rules (also included in the manuscript as Supplementary File 2).
- `retrieval_flow_diagram.png` — the sequence-retrieval flow diagram
  (Additional File 2 in the manuscript).

### `manuscript/`

- `Streptococcus_anginosus_manuscript_clean.docx` — manuscript with all
  tracked changes accepted (submission-ready main manuscript file).
- `Streptococcus_anginosus_manuscript_tracked_changes.docx` — the same
  manuscript with full tracked-changes markup against the original submission
  (for editor/reviewer reference). Page/line numbers cited in the response
  letter refer to this file, using its own built-in line numbering (restarts
  at 1 on each page).
- `Streptococcus_anginosus_manuscript_clean_colored-edits.docx` — the same
  final text as the clean file, but with every added/edited passage shown in
  blue (no strikethrough, no change bars) so a reader can see at a glance
  what changed without wading through full track-changes markup.
- `Response_to_Editor_and_Reviewers.docx` — the point-by-point response to
  the editor's decision letter and all four reviewers, with a manuscript
  page/line citation under each resolved item.

## Reproducing the analysis

```bash
git clone [AUTHOR TO CONFIRM: repository URL]
cd [AUTHOR TO CONFIRM: repository name]
pip install -r requirements.txt
cd code
python parsimony_permutation_test.py                       # full dataset
python parsimony_permutation_test.py --exclude-composition-failing
python parsimony_permutation_test.py --leave-major-study-out
python parsimony_permutation_test.py --directly-reported-only
```

Each run loads `../data/curated_metadata.csv` and
`../phylogenetics/final_tree/Aligned.Strepto.aln.treefile`, matches tree tips
to metadata by accession number, computes the observed Fitch parsimony score
for the "country" character, and compares it to a null distribution built
from 9,999 label-shuffling permutations (fixed seed 20260827 by default, for
exact reproducibility — pass `--seed` to change it).

## Requirements

Python 3.9+, `pandas`, `dendropy`, `openpyxl` (see `requirements.txt`).

## License

Code (`code/`) is released under the MIT License (see [`LICENSE`](LICENSE)).
The manuscript and curated data are not covered by that license — see the
note at the end of `LICENSE`.
