#!/usr/bin/env python3
"""
make_metadata_figures.py

Generates Figures 1-5 (metadata completeness / breakdown charts) as real
PNG images, using numbers independently re-derived from the curated
350-record dataset (curated350.csv) and cross-checked against the exact
figures already reported in the manuscript text (Results: "Metadata
completeness" and "Reported geographic origin, host, and sample source"
sections).

Verification performed before writing this script (see chat transcript):
  - Country breakdown (Figure 3) reproduced exactly via
    github_repo/code/country_mapping.py: 274 known-country records,
    19-country distribution matches the manuscript's Figure 3 description
    number-for-number.
  - Host breakdown (Figure 4) reproduced exactly by hand-classifying every
    non-missing "Host" free-text value: 289 human-associated, 1 dog
    (Canis familiaris), 1 sea water, 1 domestic wastewater = 292, matching
    the manuscript's "289 (99.0%) human-associated ... one domestic dog
    sample, one seawater sample, and one domestic wastewater sample".
  - Sample-source completeness (n=242 reporting, 108 missing, treating
    "---" and "Unknown" as missing) reproduced exactly.
  - Sample-source category breakdown (Figure 5) uses the manuscript's
    already-stated 9-category counts (125/35/34/22/11/6/5/3/1 = 242);
    two of the nine categories (Oral/Dental=125, Throat=35) and the Blood
    category (34) were independently reconstructed from the raw free-text
    values and matched exactly, and all nine sum to exactly 242 -- the
    original per-record assignment for the remaining categories was not
    re-derivable from a saved script, so those counts are taken as-is from
    the already-verified manuscript text rather than re-guessed.

No numbers are invented here; every count either comes from a script run
against the real data in this repo, or from manuscript text whose totals
were cross-checked against the real data.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

GREEN = "#2E7D32"
RED = "#C62828"
BLUE = "#1565C0"

plt.rcParams.update({
    "font.size": 11,
    "figure.facecolor": "white",
    "savefig.facecolor": "white",
})

# ---------------------------------------------------------------------
# Figure 1: bar chart, reported vs not reported, 3 fields (n=350)
# ---------------------------------------------------------------------
def figure1():
    categories = ["Sample Source", "Geographic Origin", "Host"]
    reported = [242, 274, 292]
    not_reported = [108, 76, 58]

    x = range(len(categories))
    width = 0.35
    fig, ax = plt.subplots(figsize=(9, 5.5))
    b1 = ax.bar([i - width/2 for i in x], reported, width, label="Reported", color=GREEN)
    b2 = ax.bar([i + width/2 for i in x], not_reported, width, label="Not reported", color=RED)

    for bars in (b1, b2):
        for bar in bars:
            h = bar.get_height()
            ax.annotate(f"{h}", xy=(bar.get_x() + bar.get_width()/2, h),
                        xytext=(0, 3), textcoords="offset points",
                        ha="center", va="bottom", fontsize=10)

    ax.set_xticks(list(x))
    ax.set_xticklabels(categories)
    ax.set_ylabel("Number of records")
    ax.set_title("Metadata Completeness — NCBI 16S Streptococcus anginosus Records (n=350)", fontsize=12)
    ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.02))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, 345)
    fig.tight_layout()
    fig.savefig("Figure1_metadata_completeness_bar.png", dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------
# Figure 2: three side-by-side pie charts (n=350)
# ---------------------------------------------------------------------
def figure2():
    panels = [
        ("Sample Source", 242, 108),
        ("Geographic Origin", 274, 76),
        ("Host", 292, 58),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(11, 4.2))
    for ax, (title, rep, miss) in zip(axes, panels):
        vals = [rep, miss]
        pct_rep = 100 * rep / (rep + miss)
        pct_miss = 100 - pct_rep
        wedges, texts, autotexts = ax.pie(
            vals, colors=[GREEN, RED], startangle=90, counterclock=False,
            autopct=lambda p: f"{p:.0f}%",
            wedgeprops={"edgecolor": "white", "linewidth": 1.2},
        )
        for t in autotexts:
            t.set_color("white")
            t.set_fontweight("bold")
        ax.set_title(f"{title}\n({rep} Reported, {miss} Missing)", fontsize=10)
        ax.set_aspect("equal")

    fig.legend(handles=[
                   plt.Rectangle((0, 0), 1, 1, color=GREEN),
                   plt.Rectangle((0, 0), 1, 1, color=RED),
               ], labels=["Reported", "Missing"], loc="lower center", ncol=2,
               frameon=False)
    fig.suptitle("Proportion of Records with Reported vs. Missing Metadata (n=350)", fontsize=12)
    fig.tight_layout(rect=[0, 0.06, 1, 0.93])
    fig.savefig("Figure2_metadata_completeness_pies.png", dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------
# Figure 3: horizontal bar chart, geographic origin by country (n=274)
# ---------------------------------------------------------------------
def figure3():
    data = [
        ("South Korea", 83), ("Denmark", 46), ("Spain", 37), ("China", 34),
        ("Taiwan", 17), ("Japan", 11), ("India", 10), ("Netherlands", 8),
        ("USA", 8), ("Sweden", 8), ("Poland", 3), ("Turkey", 2),
        ("Saudi Arabia", 1), ("Mexico", 1), ("UK", 1), ("Iraq", 1),
        ("Libya", 1), ("Germany", 1), ("Egypt", 1),
    ]
    assert sum(v for _, v in data) == 274
    data_sorted = sorted(data, key=lambda t: t[1])  # ascending for horizontal bar (largest at top)
    labels = [d[0] for d in data_sorted]
    values = [d[1] for d in data_sorted]

    fig, ax = plt.subplots(figsize=(8, 7.5))
    bars = ax.barh(labels, values, color=BLUE)
    for bar, v in zip(bars, values):
        ax.annotate(f"{v}", xy=(bar.get_width(), bar.get_y() + bar.get_height()/2),
                    xytext=(4, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=9)
    ax.set_xlabel("Number of records")
    ax.set_title("Geographic Origin by Country (n=274 reporting)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0, 95)
    fig.tight_layout()
    fig.savefig("Figure3_geographic_origin_by_country.png", dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------
# Figure 4: host category (n=292) -- bar chart (3 of 4 categories are n=1,
# so a bar chart communicates this far better than a pie chart would)
# ---------------------------------------------------------------------
def figure4():
    data = [
        ("Human", 289),
        ("Domestic dog\n(Canis familiaris)", 1),
        ("Sea water", 1),
        ("Domestic wastewater", 1),
    ]
    assert sum(v for _, v in data) == 292
    labels = [d[0] for d in data]
    values = [d[1] for d in data]
    colors = [GREEN, "#8D6E63", "#0288D1", "#6D4C41"]

    fig, ax = plt.subplots(figsize=(7, 5))
    bars = ax.bar(labels, values, color=colors)
    ax.set_yscale("symlog", linthresh=1)
    for bar, v in zip(bars, values):
        ax.annotate(f"{v}", xy=(bar.get_x() + bar.get_width()/2, v),
                    xytext=(0, 4), textcoords="offset points",
                    ha="center", va="bottom", fontsize=10, fontweight="bold")
    ax.set_ylabel("Number of records (log scale)")
    ax.set_title("Host Category Among Records with Reported Host Metadata\n(n=292 reporting; 289/292 = 99.0% human-associated)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(0, 500)
    ax.yaxis.set_major_formatter(mticker.ScalarFormatter())
    fig.tight_layout()
    fig.savefig("Figure4_host_category.png", dpi=300)
    plt.close(fig)


# ---------------------------------------------------------------------
# Figure 5: horizontal bar chart, sample source by clinical category (n=242)
# ---------------------------------------------------------------------
def figure5():
    data = [
        ("Oral / Dental", 125),
        ("Throat", 35),
        ("Blood", 34),
        ("Deep tissue / Abscess / Respiratory", 22),
        ("Gastrointestinal", 11),
        ("Other / Unspecified", 6),
        ("Urogenital", 5),
        ("Clinical isolate (unspecified site)", 3),
        ("Environmental", 1),
    ]
    assert sum(v for _, v in data) == 242
    data_sorted = sorted(data, key=lambda t: t[1])
    labels = [d[0] for d in data_sorted]
    values = [d[1] for d in data_sorted]

    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    bars = ax.barh(labels, values, color="#6A1B9A")
    for bar, v in zip(bars, values):
        ax.annotate(f"{v}", xy=(bar.get_width(), bar.get_y() + bar.get_height()/2),
                    xytext=(4, 0), textcoords="offset points",
                    ha="left", va="center", fontsize=9)
    ax.set_xlabel("Number of sequences")
    ax.set_title("Sample Source by Clinical Category (n=242 reporting)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(0, 140)
    fig.tight_layout()
    fig.savefig("Figure5_sample_source_by_category.png", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
    figure4()
    figure5()
    print("Wrote Figure1-5 PNGs.")
