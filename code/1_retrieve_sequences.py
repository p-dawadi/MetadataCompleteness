# 1_retrieve_sequences.py

import io
import os
import time

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
from Bio import Entrez, SeqIO
import pandas as pd

# ---- REQUIRED BY NCBI ----
EMAIL = "pdawadi@go.olemiss.edu"                 # your email, as required by NCBI's usage policy
TOOL = "streptococcus_anginosus_16S_retrieval"
API_KEY = None            # optional: put your NCBI API key string here if you have one

Entrez.email = EMAIL
Entrez.tool = TOOL
if API_KEY:
    Entrez.api_key = API_KEY

ORGANISM = "Streptococcus anginosus"

# 16S-focused query (rRNA feature or biomolecule + "16S" term) -- this is
# the ACTUAL query used for the dataset already in the manuscript.
QUERY = (
    f'("{ORGANISM}"[Organism]) AND '
    f'((rrna[Feature Key]) OR (biomol_rRNA[PROP])) AND '
    f'(16S[All Fields])'
)

# Length thresholds -- applied AFTER retrieval, in this script, not as
# part of the Entrez query string itself.
MIN_LEN_ALL = 200
MIN_LEN_FULL, MAX_LEN_FULL = 1200, 1700

ESEARCH_PAGE = 5000       # IDs per esearch page
ESUMMARY_BATCH = 500      # titles per esummary batch (small XML responses; batching is fine)
REQUEST_PAUSE = 0.34      # ~3 requests/sec (NCBI's no-API-key limit)

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
CONNECT_TIMEOUT = 10
READ_TIMEOUT = 60

CACHE_FILE = "genbank_cache.gb"
FETCHED_IDS_FILE = "fetched_ids.txt"
FAILED_IDS_FILE = "failed_ids.txt"


def make_session():
    """A requests Session with connection-level retries baked in --
    handles the interrupted/incomplete downloads your network has been
    producing far better than raw urllib does."""
    session = requests.Session()
    retry = Retry(
        total=5,
        connect=5,
        read=5,
        backoff_factor=1.5,               # 1.5s, 3s, 6s, 12s, 24s
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=frozenset(["GET", "POST"]),
        raise_on_status=False,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_maxsize=4)
    session.mount("https://", adapter)
    session.mount("http://", adapter)
    return session


SESSION = make_session()


def eutils_params(extra):
    params = {"email": EMAIL, "tool": TOOL}
    if API_KEY:
        params["api_key"] = API_KEY
    params.update(extra)
    return params


def fetch_all_ids(query, page=ESEARCH_PAGE):
    ids, retstart = [], 0
    while True:
        params = eutils_params({
            "db": "nucleotide", "term": query, "retmax": page,
            "retstart": retstart, "retmode": "json",
        })
        r = SESSION.get(f"{EUTILS_BASE}/esearch.fcgi", params=params,
                         timeout=(CONNECT_TIMEOUT, READ_TIMEOUT))
        r.raise_for_status()
        chunk = r.json().get("esearchresult", {}).get("idlist", [])
        if not chunk:
            break
        ids.extend(chunk)
        if len(chunk) < page:
            break
        retstart += page
        time.sleep(REQUEST_PAUSE)
    return ids


def fetch_titles_by_accession(uids):
    """Return dict: accession(without version) -> Title, via ESummary.

    Uses POST (form body) rather than GET (query string) -- with 500 IDs
    per batch the URL gets long enough that some networks/proxies (or
    NCBI itself) reject it with "414 Request-URI Too Long". A POST body
    has no such limit and is what NCBI itself recommends for large ID
    lists."""
    acc_to_title = {}
    for i in range(0, len(uids), ESUMMARY_BATCH):
        chunk = uids[i:i + ESUMMARY_BATCH]
        params = eutils_params({
            "db": "nucleotide", "id": ",".join(chunk), "retmode": "json",
        })
        r = SESSION.post(f"{EUTILS_BASE}/esummary.fcgi", data=params,
                          timeout=(CONNECT_TIMEOUT, READ_TIMEOUT))
        r.raise_for_status()
        result = r.json().get("result", {})
        for uid in result.get("uids", []):
            doc = result.get(uid, {})
            caption = str(doc.get("caption", "")).strip()
            title = str(doc.get("title", "")).strip()
            if caption:
                acc = caption.split(".")[0]
                if title:
                    acc_to_title[acc] = title
        time.sleep(REQUEST_PAUSE)
    return acc_to_title


def load_fetched_ids():
    if not os.path.exists(FETCHED_IDS_FILE):
        return set()
    with open(FETCHED_IDS_FILE) as f:
        return {line.strip() for line in f if line.strip()}


def fetch_and_cache_all(ids):
    """Fetch each record ONE AT A TIME, appending successes immediately
    to CACHE_FILE and recording done IDs in FETCHED_IDS_FILE. Safe to
    interrupt and re-run -- already-fetched IDs are skipped."""
    already_done = load_fetched_ids()
    remaining = [uid for uid in ids if uid not in already_done]
    print(f"{len(already_done)} of {len(ids)} records already cached "
          f"(from a previous run); {len(remaining)} left to fetch.")

    failed = []
    with open(CACHE_FILE, "a") as cache_f, \
         open(FETCHED_IDS_FILE, "a") as done_f:
        for n, uid in enumerate(remaining, 1):
            ok = False
            last_exc = None
            for attempt in range(1, 4):
                try:
                    params = eutils_params({
                        "db": "nucleotide", "id": uid,
                        "rettype": "gb", "retmode": "text",
                    })
                    r = SESSION.get(f"{EUTILS_BASE}/efetch.fcgi", params=params,
                                     timeout=(CONNECT_TIMEOUT, READ_TIMEOUT))
                    r.raise_for_status()
                    text = r.text
                    if not text.strip():
                        raise ValueError("empty response body")
                    cache_f.write(text.rstrip("\n") + "\n")
                    cache_f.flush()
                    done_f.write(uid + "\n")
                    done_f.flush()
                    ok = True
                    break
                except Exception as exc:
                    last_exc = exc
                    time.sleep(1.5 * attempt)
            if not ok:
                print(f"  [FAILED] id={uid} after 3 attempts: {last_exc}")
                failed.append(uid)
            if n % 50 == 0 or n == len(remaining):
                print(f"  ... {n}/{len(remaining)} fetched this run "
                      f"({len(already_done) + n - len(failed)} total cached, "
                      f"{len(failed)} failed)")
            time.sleep(REQUEST_PAUSE)

    if failed:
        with open(FAILED_IDS_FILE, "a") as f:
            for uid in failed:
                f.write(uid + "\n")
        print(f"{len(failed)} records failed all retries this run -- logged to "
              f"{FAILED_IDS_FILE}. Just re-run the script to retry them "
              f"(everything else stays cached).")


def extract_source(rec):
    """Return isolation/source-like text from the source feature (best-effort)."""
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


def save_fasta_by_accessions(acc_list, outfile):
    """Fetch FASTA one accession at a time (same fireproof pattern) --
    guarantees concrete sequences and survives the same flaky-network
    failure mode as the GenBank fetch above."""
    acc_list = list(dict.fromkeys([a.split(".")[0] for a in acc_list]))  # de-dup, strip version
    done_file = outfile + ".done.txt"
    already_done = set()
    if os.path.exists(done_file):
        with open(done_file) as f:
            already_done = {line.strip() for line in f if line.strip()}

    with open(outfile, "a") as out, open(done_file, "a") as done_f:
        remaining = [a for a in acc_list if a not in already_done]
        print(f"FASTA: {len(already_done)} of {len(acc_list)} already saved to {outfile}; "
              f"{len(remaining)} left.")
        for n, acc in enumerate(remaining, 1):
            for attempt in range(1, 4):
                try:
                    params = eutils_params({
                        "db": "nucleotide", "id": acc,
                        "rettype": "fasta", "retmode": "text",
                    })
                    r = SESSION.get(f"{EUTILS_BASE}/efetch.fcgi", params=params,
                                     timeout=(CONNECT_TIMEOUT, READ_TIMEOUT))
                    r.raise_for_status()
                    text = r.text
                    if not text.strip():
                        raise ValueError("empty response body")
                    out.write(text.rstrip("\n") + "\n")
                    out.flush()
                    done_f.write(acc + "\n")
                    done_f.flush()
                    break
                except Exception as exc:
                    if attempt == 3:
                        print(f"  [FAILED] {acc}: {exc}")
                    time.sleep(1.5 * attempt)
            if n % 50 == 0 or n == len(remaining):
                print(f"  ... {n}/{len(remaining)} fetched this run for {outfile}")
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

    # Titles (small ESummary calls, batched -- not the failure point so
    # far, left batched for speed).
    print("Fetching titles via ESummary...")
    title_map = fetch_titles_by_accession(ids)

    # The step that was crashing: now one-record-at-a-time, cached,
    # resumable.
    print("Fetching full GenBank records (one at a time, cached, resumable)...")
    fetch_and_cache_all(ids)

    if os.path.exists(FAILED_IDS_FILE):
        with open(FAILED_IDS_FILE) as f:
            n_failed = len({line.strip() for line in f if line.strip()})
        if n_failed:
            print(f"NOTE: {n_failed} record(s) are still unfetched after retries "
                  f"(see {FAILED_IDS_FILE}). Re-run this script to retry just those; "
                  f"the summary below only reflects what's currently cached.")

    print("Parsing cached GenBank records...")
    rows_all, rows_full = [], []
    with open(CACHE_FILE) as f:
        cache_text = f.read()
    for rec in SeqIO.parse(io.StringIO(cache_text), "genbank"):
        acc = (rec.annotations.get("accessions", [rec.id])[0] or rec.id).split(".")[0]
        title = title_map.get(acc) or extract_definition(rec)
        length = safe_len(rec)
        source = extract_source(rec)
        ncbi_link = f"https://www.ncbi.nlm.nih.gov/nuccore/{acc}"

        row = {
            "Accession Number": acc,
            "Title": title,
            "Sequence Length": length,
            "Source of sample": source,
            "NCBI Link": ncbi_link,
        }

        if length >= MIN_LEN_ALL:
            rows_all.append(row)
        if MIN_LEN_FULL <= length <= MAX_LEN_FULL:
            rows_full.append(row)

    csv_all = "Sp_all_ge200bp.csv"
    csv_full = "Sp_1200_1700bp.csv"
    fa_all = "Sp_16s_ge200bp.fasta"
    fa_full = "Sp_16S_full_1200_1700bp.fasta"

    df_all = finalize(rows_all, csv_all)
    df_full = finalize(rows_full, csv_full)

    if not df_all.empty:
        save_fasta_by_accessions(df_all["Accession Number"].tolist(), fa_all)
    if not df_full.empty:
        save_fasta_by_accessions(df_full["Accession Number"].tolist(), fa_full)

    print(f"Saved {len(df_all)} sequences to {csv_all} and {fa_all}.")
    print(f"Saved {len(df_full)} sequences to {csv_full} and {fa_full}.")


if __name__ == "__main__":
    main()
