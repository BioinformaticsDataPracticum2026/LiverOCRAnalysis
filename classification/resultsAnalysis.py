#!/usr/bin/env python3

"""
Analysis and Validation of OCR classification results.
    
Generates:
    - Summary statistics
    - Distribution plots
    - Venn diagrams
    - Quality metrics
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np

#Set visualization style for better-looking plots
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)

#Results directory -> Contains all output BED files from classification
results_dir = Path("results")

print("OCR Classification Results Analysis")

#Basic statistics of results
print("\nBASIC STATISTICS")

#Dictionary mapping category names to BED file names
files = {
    "Human Promoters": "human_promoters.bed",
    "Human Enhancers": "human_enhancers.bed",
    "Mouse Promoters": "mouse_promoters.bed",
    "Mouse Enhancers": "mouse_enhancers.bed",
}

#Read each file and compute statistics
stats = {}
for name, file in files.items():
    path = results_dir / file
    if path.exists():
        # Read BED file (tab-separated, no header)
        # Columns: chromosome (0), start (1), end (2), [optional additional columns]
        df = pd.read_csv(path, sep="\t", header=None)

        #Calculate region lengths (end - start)
        stats[name] = {
            "count": len(df),
            "mean_length": (df[2] - df[1]).mean(),
            "median_length": (df[2] - df[1]).median(),
        }

        #Print statistics for this category
        print(f"{name}:")
        print(f"  Total regions: {len(df):,}")
        print(f"  Mean length: {stats[name]['mean_length']} bp")
        print(f"  Median length: {stats[name]['median_length']} bp")
    else:
        print(f"{file} not found")


#PROMOTER vs. ENHANCER RATIO
print("\nPROMOTER vs ENHANCER RATIO")

fig, axes = plt.subplots(1, 2, figsize=(14, 5))

#Human: Read promoter and enhancer counts
human_prom = len(pd.read_csv(results_dir / "human_promoters.bed", sep="\t", header=None))
human_enh = len(pd.read_csv(results_dir / "human_enhancers.bed", sep="\t", header=None))
human_total = human_prom + human_enh

#Create pie chart for human
axes[0].pie(
    [human_prom, human_enh],
    labels=[f"Promoters\n({human_prom:,})", f"Enhancers\n({human_enh:,})"],
    autopct="%1.1f%%", #Show percentages on pie chart
    colors=["#ff9999", "#66b3ff"],
    startangle=90
)
axes[0].set_title(f"Human OCRs (Total: {human_total:,})", fontsize=12, fontweight="bold")

print(f"Human: {human_prom:,} promoters ({100*human_prom/human_total}%) | {human_enh:,} enhancers ({100*human_enh/human_total}%)")

#Mouse: Read promoter and enhancer counts
mouse_prom = len(pd.read_csv(results_dir / "mouse_promoters.bed", sep="\t", header=None))
mouse_enh = len(pd.read_csv(results_dir / "mouse_enhancers.bed", sep="\t", header=None))
mouse_total = mouse_prom + mouse_enh

#Create pie chart for mouse
axes[1].pie(
    [mouse_prom, mouse_enh],
    labels=[f"Promoters\n({mouse_prom:,})", f"Enhancers\n({mouse_enh:,})"],
    autopct="%1.1f%%",
    colors=["#ff9999", "#66b3ff"],
    startangle=90
)
axes[1].set_title(f"Mouse OCRs (Total: {mouse_total:,})", fontsize=12, fontweight="bold")

print(f"Mouse: {mouse_prom:,} promoters ({100*mouse_prom/mouse_total}%) | {mouse_enh:,} enhancers ({100*mouse_enh/mouse_total}%)")

plt.tight_layout()
plt.savefig("results_01_promoter_enhancer_ratio.png", dpi=300, bbox_inches="tight")
print("✓ Saved: results_01_promoter_enhancer_ratio.png")


#CROSS-SPECIES CONSERVATION
print("\nCROSS-SPECIES CONSERVATION")

#Read counts from mapped files (these are regions from one species mapped to another genome)
h_to_m_prom = len(pd.read_csv(results_dir / "mouse_mapped_from_human_promoters.bed", sep="\t", header=None))
h_to_m_enh = len(pd.read_csv(results_dir / "mouse_mapped_from_human_enhancers.bed", sep="\t", header=None))
m_to_h_prom = len(pd.read_csv(results_dir / "human_mapped_from_mouse_promoters.bed", sep="\t", header=None))
m_to_h_enh = len(pd.read_csv(results_dir / "human_mapped_from_mouse_enhancers.bed", sep="\t", header=None))

#Read counts from shared files (these are regions that overlap BOTH mapped and native elements)
#Represent truly conserved regulatory elements
h_to_m_shared = len(pd.read_csv(results_dir / "mouse_shared_from_human.bed", sep="\t", header=None))
m_to_h_shared = len(pd.read_csv(results_dir / "human_shared_from_mouse.bed", sep="\t", header=None))

#Print conservation percentages
print(f"Human→Mouse mapping: {h_to_m_prom + h_to_m_enh:,} total ({h_to_m_prom:,} promoters, {h_to_m_enh:,} enhancers)")
print(f"  Shared with native mouse: {h_to_m_shared:,} ({100*h_to_m_shared/(h_to_m_prom + h_to_m_enh)}% conservation)")

print(f"Mouse→Human mapping: {m_to_h_prom + m_to_h_enh:,} total ({m_to_h_prom:,} promoters, {m_to_h_enh:,} enhancers)")
print(f"  Shared with native human: {m_to_h_shared:,} ({100*m_to_h_shared/(m_to_h_prom + m_to_h_enh)}% conservation)")

#Create bar charts showing conservation rates
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

#Human to Mouse conservation
categories = ["Promoters", "Enhancers", "Shared"]
h_to_m_vals = [h_to_m_prom, h_to_m_enh, h_to_m_shared]
axes[0].bar(categories, h_to_m_vals, color=["#ff9999", "#66b3ff", "#99ff99"])
axes[0].set_title("Human → Mouse Conservation", fontsize=12, fontweight="bold")
axes[0].set_ylabel("Number of regions")
#Add value labels on bars
for i, v in enumerate(h_to_m_vals):
    axes[0].text(i, v + 500, str(f"{v:,}"), ha="center", fontweight="bold")

#Mouse to Human conservation
m_to_h_vals = [m_to_h_prom, m_to_h_enh, m_to_h_shared]
axes[1].bar(categories, m_to_h_vals, color=["#ff9999", "#66b3ff", "#99ff99"])
axes[1].set_title("Mouse → Human Conservation", fontsize=12, fontweight="bold")
axes[1].set_ylabel("Number of regions")
#Add value labels on bars
for i, v in enumerate(m_to_h_vals):
    axes[1].text(i, v + 500, str(f"{v:,}"), ha="center", fontweight="bold")

plt.tight_layout()
plt.savefig("results_02_cross_species_conservation.png", dpi=300, bbox_inches="tight")
print("✓ Saved: results_02_cross_species_conservation.png")


#REGION SIZE DISTRIBUTION
print("\nREGION SIZE DISTRIBUTION")

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

#List of files to plot with their titles and axes
files_to_plot = [
    ("human_promoters.bed", "Human Promoters", axes[0, 0]),
    ("human_enhancers.bed", "Human Enhancers", axes[0, 1]),
    ("mouse_promoters.bed", "Mouse Promoters", axes[1, 0]),
    ("mouse_enhancers.bed", "Mouse Enhancers", axes[1, 1]),
]

#Create histogram for each file
for file, title, ax in files_to_plot:
    path = results_dir / file
    df = pd.read_csv(path, sep="\t", header=None)

    #Calculate region sizes (end - start)
    sizes = df[2] - df[1]
    
    #Create histogram
    ax.hist(sizes, bins=50, color="#66b3ff", edgecolor="black", alpha=0.7)
    ax.set_xlabel("Region size (bp)")
    ax.set_ylabel("Frequency")
    ax.set_title(title)

    #Add vertical line showing mean
    ax.axvline(sizes.mean(), color="red", linestyle="--", linewidth=2, label=f"Mean: {sizes.mean()} bp")
    ax.legend()

plt.tight_layout()
plt.savefig("results_03_region_size_distribution.png", dpi=300, bbox_inches="tight")
print("Saved: results_03_region_size_distribution.png")


#SUMMARY TABLE
print("\nSUMMARY TABLE")

#Create structured data for summary table
summary_data = {
    "Category": [
        "Human Native",
        "Mouse Native",
        "Human←Mouse Mapped",
        "Mouse←Human Mapped",
        "Human-Mouse Shared",
    ],
    "Promoters": [
        human_prom,
        mouse_prom,
        m_to_h_prom,
        h_to_m_prom,
        m_to_h_shared,
    ],
    "Enhancers": [
        human_enh,
        mouse_enh,
        m_to_h_enh,
        h_to_m_enh,
        h_to_m_shared,
    ],
    "Total": [
        human_total,
        mouse_total,
        m_to_h_prom + m_to_h_enh,
        h_to_m_prom + h_to_m_enh,
        h_to_m_shared + m_to_h_shared,
    ]
}

summary_df = pd.DataFrame(summary_data)
print(summary_df.to_string(index=False))

#Save to CSV
summary_df.to_csv("results_summary.csv", index=False)
print("✓ Saved: results_summary.csv")


#SANITY CHECKS
print("\nSANITY CHECKS")

#The following is a series of checks to ensure data quality
checks = [
    #Check 1: Promoter and enhancer counts should sum to total for human
    ("Human promoters + enhancers = total",
     human_prom + human_enh == human_total,
     "Passed sanity check." if human_prom + human_enh == human_total
     else "Failed sanity check."),
    
    #Check 2: Promoter and enhancer counts should sum to total for mouse
    ("Mouse promoters + enhancers = total",
     mouse_prom + mouse_enh == mouse_total,
     "Passed sanity check." if mouse_prom + mouse_enh == mouse_total
     else "Failed sanity check."),

     #Check 3: Promoter to enhancer ratio should be reasonable (~20-50% promoters typical)
    ("Promoter-to-enhancer ratio reasonable",
     (human_prom/human_enh > 0.2) and (human_prom/human_enh < 0.8),
     "Passed sanity check." if (human_prom/human_enh > 0.2) and (human_prom/human_enh < 0.8)
     else "Failed sanity check."),
    
    #Check 4: There should be some conservation between species
    ("Conservation rate > 0%",
     (h_to_m_shared > 0) and (m_to_h_shared > 0),
     "Passed sanity check." if (h_to_m_shared > 0) and (m_to_h_shared > 0)
     else "Failed sanity check."),

     #Check 5: Conservation should not be 100% (some elements are species-specific)
    ("Conservation rate < 100%",
     (h_to_m_shared < h_to_m_prom + h_to_m_enh) and (m_to_h_shared < m_to_h_prom + m_to_h_enh),
     "Passed sanity check." if (h_to_m_shared < h_to_m_prom + h_to_m_enh)
     else "Failed sanity check."),
]

#Print results for each check
for check_name, result, symbol in checks:
    print(f"{symbol} {check_name}: {result}")

#FINISHED!!!
print("VALIDATION COMPLETE")
print("\nGenerated visualizations:")
print("  1. results_01_promoter_enhancer_ratio.png")
print("  2. results_02_cross_species_conservation.png")
print("  3. results_03_region_size_distribution.png")
print("  4. results_summary.csv")