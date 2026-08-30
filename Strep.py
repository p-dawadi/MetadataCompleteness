# lgarvieae_16S_fetch_mincols.py
# pip install biopython pandas
from Bio import Entrez, SeqIO
import pandas as pd
import time

# ---- REQUIRED BY NCBI ----
Entrez.email = "pdawadi@go.olemiss.edu"          # <-- put your email
Entrez.tool  = "lactococcus_garvieae_downloader"
# Optional: Entrez.api_key = "YOUR_NCBI_API_KEY"

ORGANISM = "Streptococcus anginosus"

# 16S-focused query (rRNA feature or biomolecule + "16S" term)
QUERY = (
    f'("{ORGANISM}"[Organism]) AND '
    f'((rrna[Feature Key]) OR (biomol_rRNA[PROP])) AND '
    f'(16S[All Fields])'
)

# length thresholds
MIN_LEN_ALL = 200
MIN_LEN_FULL, MAX_LEN_FULL = 1200, 1700

# API hygiene
ESEARCH_PAGE  = 5000     # IDs per esearch page
EFETCH_BATCH  = 300      # records per efetch batch (<=500 is fine)
ESUMMARY_BATCH = 500     # titles per esummary batch
REQUEST_PAUSE = 0.34     # ~3 requests/sec

# ---------- helpers ----------
def fetch_all_ids(query, page=ESEARCH_PAGE):
    ids, retstart = [], 0
    while True:
        h = Entrez.esearch(db="nucleotide", term=query, retmax=page, retstart=retstart)
        rec = Entrez.read(h); h.close()
        chunk = rec.get("IdList", [])
        if not chunk:
            break
        ids.extend(chunk)
        if len(chunk) < page:
            break
        retstart += page
        time.sleep(REQUEST_PAUSE)
    return ids

def fetch_titles_by_accession(uids):
    """
    Return dict: accession(without version) -> Title
    Uses ESummary fields: Caption (accession.version) and Title.
    """
    acc_to_title = {}
    for i in range(0, len(uids), ESUMMARY_BATCH):
        chunk = uids[i:i+ESUMMARY_BATCH]
        h = Entrez.esummary(db="nucleotide", id=",".join(chunk), retmode="xml")
        summ = Entrez.read(h); h.close()
        for doc in summ:
            # Caption is accession.version; strip version to match our 'acc'
            caption = str(doc.get("Caption", "")).strip()
            title   = str(doc.get("Title", "")).strip()
            if caption:
                acc = caption.split(".")[0]
                if title:
                    acc_to_title[acc] = title
        time.sleep(REQUEST_PAUSE)
    return acc_to_title

def extract_source(rec):
    """Return isolation/source-like text from source feature (best-effort)."""
    for feat in rec.features:
        if feat.type == "source":
            q = feat.qualifiers
            return (
                q.get("isolation_source", [""]) or
                q.get("tissue_type", [""]) or
                q.get("host", [""])
            )[0]
    return ""

def extract_definition(rec):
    """Fallback: GenBank DEFINITION line parsed by Biopython."""
    try:
        return (rec.description or "").strip()
    except Exception:
        return ""

def safe_len(rec):
    """Robust length even if sequence is UnknownSeq."""
    try:
        if rec.seq and rec.seq.__class__.__name__ != "UnknownSeq":
            return len(rec.seq)
    except Exception:
        pass
    ann_len = rec.annotations.get("sequence_length")
    try:
        return int(ann_len)
    except Exception:
        return 0

def save_fasta_by_accessions(acc_list, outfile, batch=EFETCH_BATCH):
    """Fetch FASTA by accession strings (guarantees concrete sequences)."""
    acc_list = list(dict.fromkeys([a.split(".")[0] for a in acc_list]))  # de-dup, strip version
    with open(outfile, "w") as out:
        for i in range(0, len(acc_list), batch):
            chunk = acc_list[i:i+batch]
            h = Entrez.efetch(db="nucleotide", id=",".join(chunk),
                              rettype="fasta", retmode="text")
            out.write(h.read().strip() + "\n")
            h.close()
            time.sleep(REQUEST_PAUSE)

def finalize(rows, csv_name):
    cols = ["Accession Number", "Title", "Sequence Length", "Source of sample", "NCBI Link"]
    df = pd.DataFrame(rows, columns=cols)
    if df.empty:
        df.to_csv(csv_name, index=False)
        return df
    df = (df
          .drop_duplicates(subset=["Accession Number"])
          .sort_values("Sequence Length", ascending=False)
          .reset_index(drop=True))
    df.to_csv(csv_name, index=False)
    return df

# ---------- main ----------
def main():
    ids = fetch_all_ids(QUERY)
    print(f"Found {len(ids)} nucleotide IDs for 16S rRNA of {ORGANISM}")

    rows_all, rows_full = [], []

    # Process in the same EFETCH batches; get titles for each batch via ESummary
    for i in range(0, len(ids), EFETCH_BATCH):
        uid_batch = ids[i:i+EFETCH_BATCH]

        # 1) Titles by accession (from ESummary)
        title_map = fetch_titles_by_accession(uid_batch)  # acc -> Title

        # 2) Full GenBank records
        h = Entrez.efetch(db="nucleotide", id=uid_batch, rettype="gb", retmode="text")
        for rec in SeqIO.parse(h, "genbank"):
            acc = (rec.annotations.get("accessions", [rec.id])[0] or rec.id).split(".")[0]
            title = title_map.get(acc) or extract_definition(rec)  # prefer ESummary Title
            length = safe_len(rec)
            source = extract_source(rec)
            ncbi_link = f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}"

            row = {
                "Accession Number": acc,
                "Title": title,
                "Sequence Length": length,
                "Source of sample": source,
                "NCBI Link": ncbi_link
            }

            if length >= MIN_LEN_ALL:
                rows_all.append(row)
            if MIN_LEN_FULL <= length <= MAX_LEN_FULL:
                rows_full.append(row)
        h.close()
        time.sleep(REQUEST_PAUSE)

    # Windows-safe filenames
    csv_all   = "Sp_all_ge200bp.csv"
    csv_full  = "Sp_1200_1700bp.csv"
    fa_all    = "Sp_16s_ge200bp.fasta"
    fa_full   = "Sp_16S_full_1200_1700bp.fasta"

    df_all  = finalize(rows_all,  csv_all)
    df_full = finalize(rows_full, csv_full)

    # FASTA pass (separate; avoids UnknownSeq issues)
    if not df_all.empty:
        save_fasta_by_accessions(df_all["Accession Number"].tolist(),  fa_all)
    if not df_full.empty:
        save_fasta_by_accessions(df_full["Accession Number"].tolist(), fa_full)

    print(f"Saved {len(df_all)} sequences to {csv_all} and {fa_all}.")
    print(f"Saved {len(df_full)} sequences to {csv_full} and {fa_full}.")

if __name__ == "__main__":
    main()
