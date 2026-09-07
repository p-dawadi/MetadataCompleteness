"""
country_mapping.py

Maps the free-text "Geographical location" field in the curated metadata
(data/curated_metadata.csv) to a standardized country name.

Two kinds of entries exist:
  1. A country stated directly ("South Korea", "China: Shenzhen").
  2. A submitting institution/hospital/university name with no country
     stated directly ("National Taiwan University Hospital"). These are
     mapped to a country by keyword match and flagged by
     `is_institution_inferred` so they can be excluded for the
     directly-reported-only sensitivity analysis (see
     parsimony_permutation_test.py --directly-reported-only).

Placeholder/missing values ("---", "unknown", "confusioninpaper", etc.)
map to None.
"""
import pandas as pd

MISSING_VALUES = {
    "---", "-", "", "-----", "unknown", "Unknown",
    "confusioninpaper", "Confusion", "Confusion ",
}

# (keyword to match, case-insensitive) -> standardized country name.
# Order matters: more specific keys should come first if they overlap.
COUNTRY_KEYWORDS = [
    ("south korea", "South Korea"), ("korea", "South Korea"),
    ("denmark", "Denmark"), ("spain", "Spain"), ("china", "China"),
    ("taiwan", "Taiwan"), ("india", "India"), ("netherlands", "Netherlands"),
    ("sweden", "Sweden"), ("usa", "USA"), ("united states", "USA"),
    ("japan", "Japan"), ("turkey", "Turkey"), ("poland", "Poland"),
    ("uk", "UK"), ("united kingdom", "UK"), ("saudi arabia", "Saudi Arabia"),
    ("germany", "Germany"), ("libya", "Libya"), ("mexico", "Mexico"),
    ("egypt", "Egypt"), ("babylon", "Iraq"),
]

# A raw value equal to one of these (verbatim) is a directly-reported
# country name, never an institution inference.
DIRECT_COUNTRY_NAMES = {
    "South Korea", "Korea", "Denmark", "Spain", "China", "India",
    "Netherlands", "Sweden", "USA", "Japan", "UK", "Saudi Arabia",
    "Libya", "Poland", "Mexico", "Egypt",
}

INSTITUTION_KEYWORDS = [
    "University", "Hospital", "Institute", "College", "Committee",
    "Sample:", "dentistry",
]


def map_country(raw):
    """Return the standardized country name for a raw metadata value, or
    None if missing/unusable."""
    if pd.isna(raw):
        return None
    s = str(raw).strip()
    if s in MISSING_VALUES:
        return None
    low = s.lower()
    for keyword, country in COUNTRY_KEYWORDS:
        if keyword in low:
            return country
    return None


def is_institution_inferred(raw):
    """True if the raw value names an institution rather than stating a
    country directly (used to build the directly-reported-only subset)."""
    if pd.isna(raw):
        return False
    s = str(raw).strip()
    if s in DIRECT_COUNTRY_NAMES:
        return False
    return any(k in s for k in INSTITUTION_KEYWORDS)


def load_curated_with_country(csv_path="../data/curated_metadata.csv"):
    """Load the curated metadata and add `country` and
    `is_institution_inferred` columns."""
    df = pd.read_csv(csv_path)
    df["country"] = df["Geographical location"].apply(map_country).astype(object)
    df["country"] = df["country"].where(pd.notna(df["country"]), None)
    df["is_institution_inferred"] = df["Geographical location"].apply(is_institution_inferred)
    return df


if __name__ == "__main__":
    df = load_curated_with_country()
    known = df[df["country"].notna()]
    print(f"Known-country records: {len(known)} of {len(df)}")
    print(known["country"].value_counts().to_string())
    inst = known[known["is_institution_inferred"]]
    print(f"\nInstitution-inferred: {len(inst)} of {len(known)} known-country records")
    print(inst["Geographical location"].value_counts().to_string())
