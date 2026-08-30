# Data dictionary

Describes [`data/curated_metadata.csv`](../data/curated_metadata.csv) /
[`data/curated_metadata.xlsx`](../data/curated_metadata.xlsx) (350 curated
*Streptococcus anginosus* 16S rRNA gene records, 1200–1700 bp, retrieved from
NCBI Nucleotide) and the derived fields produced by
[`code/country_mapping.py`](../code/country_mapping.py).

## Source columns (as curated)

| Column | Description |
|---|---|
| `Accession Number` | NCBI GenBank/RefSeq accession (version suffix, e.g. `.1`, stripped for tree-tip matching — see `base_accession()` in `code/parsimony_permutation_test.py`). |
| `Title` | Full NCBI record title as returned by the Nucleotide database (organism, strain/isolate, gene, description). |
| `Sequence Length` | Length of the deposited sequence in bp (inclusion window: 1200–1700 bp). |
| `Source of sample` | Free-text anatomical/clinical source as stated in the original GenBank record or source publication (e.g. "throat", "oral cavity", "Subgingival plaque", "Blood"); not standardized/controlled-vocabulary, and `---` where not stated. |
| `NCBI Link` | Direct URL to the NCBI Nucleotide record. |
| `Geographical location` | Free-text geographic origin as stated in the GenBank record or source publication — either a country/region stated directly, or (in many records) only a submitting institution name. See "Derived: `country`" below for how this is standardized. |
| `Host` | Free-text host/subject description as stated (e.g. "Human", "Patient", "IBD patients"); not standardized, and `---` where not stated. |
| `Article name` | Short name/identifier for the source publication or submitting study that deposited the record; used to group records by study of origin (see leave-major-study-out sensitivity analysis). |

## Missing-value convention

Placeholder strings `---`, `-`, `-----`, `unknown`/`Unknown`, `confusioninpaper`,
`Confusion` (and blank) all represent "not stated / not usable" and are treated
as missing (mapped to `None`/`NaN`), never as a category of their own.

## Derived fields (added by `code/country_mapping.py`, not present in the raw file)

| Column | Description |
|---|---|
| `country` | Standardized country name derived from `Geographical location` by case-insensitive keyword match (see `COUNTRY_KEYWORDS` in `code/country_mapping.py`). `None` when `Geographical location` is missing/a placeholder or does not match any known keyword. Used as the categorical character (the "state") in the Fitch parsimony/permutation test. |
| `is_institution_inferred` | `True` when `country` was inferred from a submitting institution/hospital/university name (matched via `INSTITUTION_KEYWORDS`: "University", "Hospital", "Institute", "College", "Committee", "Sample:", "dentistry") rather than a country stated directly in `Geographical location`. Used to build the directly-reported-only sensitivity subset (`--directly-reported-only`), which excludes all institution-inferred records and retains only records where a country/region was stated directly. |

## Country-keyword mapping table

The exact keyword → country mapping used (case-insensitive substring match on
`Geographical location`; order matters — see `COUNTRY_KEYWORDS` in
`code/country_mapping.py` for the authoritative, ordered list):

South Korea/Korea → South Korea; Denmark → Denmark; Spain → Spain; China →
China; Taiwan → Taiwan; India → India; Netherlands → Netherlands; Sweden →
Sweden; USA/United States → USA; Japan → Japan; Turkey → Turkey; Poland →
Poland; UK/United Kingdom → UK; Saudi Arabia → Saudi Arabia; Germany →
Germany; Libya → Libya; Mexico → Mexico; Egypt → Egypt; Babylon → Iraq.

A raw value that is *verbatim equal* to one of `South Korea, Korea, Denmark,
Spain, China, India, Netherlands, Sweden, USA, Japan, UK, Saudi Arabia,
Libya, Poland, Mexico, Egypt` is always treated as directly reported
(`is_institution_inferred = False`), even though it also matches a keyword.

## Related files in this folder

- `institution_to_country_mapping.csv` — the 37 raw `Geographical location`
  strings that name a submitting institution (not a country directly),
  mapped to a standardized country, with per-string record counts.
- `Supplementary_File_2_Data_Dictionary.xlsx` — the manuscript-facing
  version of this data dictionary and the institution-to-country mapping
  table, plus the explicit sample-source category rules and the three
  individually identified non-human/environmental records (one canine
  host, accession JN713285; two environmental sources — a domestic
  wastewater record and a seawater record).
- `retrieval_flow_diagram.png` — sequence-retrieval flow diagram
  (Additional File 2 in the manuscript): NCBI query → 350 records
  identified → 0 excluded on eligibility grounds → 350 included.

## Metadata-curation process notes

Conflicts between GenBank and publication metadata, when they arose, were
resolved by majority vote among the co-authors. Curation was independently
checked by four co-authors (PD, SN, BP, and MR).

## Species-identification check (type-strain alignment comparison)

The curated dataset happens to include the *S. anginosus* type strain
itself, under two accessions already present in the final alignment:
`AF104678.1` (ATCC 33397) and `NR_041722.2` (SK52 = DSM 20563, an NCBI
RefSeq targeted-locus record). `code/type_strain_identity_check.py`
computes the ungapped pairwise percent identity of every other curated
sequence against each of these directly from the alignment, as a real,
network-free substitute for BLAST/type-strain verification. Result: mean
98.4% identity (median 99.4%); 8 of 353 records fall below 95% identity,
of which the 4 lowest (42.3–51.8%) are exactly the four WGS-contig,
composition-test-failing sequences described above. See
[`results/type_strain_identity_check.txt`](../results/type_strain_identity_check.txt).

## Known limitations

- `Source of sample` and `Host` are free text as originally deposited/reported
  and have not been mapped to a controlled vocabulary in this dataset.
- Four accessions (`NZ_JAQMKN010000038`, `QSDO01000081`, `VYVX01000034`,
  `PVSY01000043`) are whole-genome-shotgun assembly contigs (SPAdes/scaffold-
  derived; see their `Title` field) rather than direct 16S rRNA gene
  amplicon/Sanger sequences; they failed IQ-TREE's composition-homogeneity
  test and are excluded in the `--exclude-composition-failing` sensitivity
  analysis (see [`results/parsimony_test_results.md`](../results/parsimony_test_results.md)).
- Country assignment reflects the geographic location *stated in the GenBank
  record or source publication*, which for many records is the submitting
  institution's location rather than an independently confirmed patient
  origin; this is exactly what `is_institution_inferred` flags, and the
  directly-reported-only sensitivity analysis addresses it directly.
