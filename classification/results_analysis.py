#!/usr/bin/env python3

"""
Analysis and Validation of OCR classification results.

Generates Output Files:
    - ortholog_status.png - Pie charts showing open vs closed orthologs
    - promoter_enhancer_classification.png - Pie charts of promoter/enhancer ratios
    - conservation_comparison.png - Stacked bar chart of conservation status
    - percentage_comparison.png - Side-by-side percentage comparison bars
    - summary_table.csv - Summary statistics table
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
results_dir = Path("raw_results")

print("OCR Classification Results Analysis")

#IDENTIFY OCRs WITH OPEN VS CLOSED ORTHOLOGS
print("\nIdentify OCRs with Open vs Closed Orthologs")

#Count regions by ortholog status
human_shared = len(pd.read_csv(results_dir / "human_shared.bed", sep="\t", header=None))
human_specific = len(pd.read_csv(results_dir / "human_specific.bed", sep="\t", header=None))
mouse_shared = len(pd.read_csv(results_dir / "mouse_shared.bed", sep="\t", header=None))
mouse_specific = len(pd.read_csv(results_dir / "mouse_specific.bed", sep="\t", header=None))

print("\nHuman OCRs:")
print(f"  With OPEN mouse orthologs:   {human_shared:,} regions")
print(f"  With CLOSED mouse orthologs: {human_specific:,} regions")
print(f"  Total classified:            {human_shared + human_specific:,} regions")

print("\nMouse OCRs:")
print(f"  With OPEN human orthologs:   {mouse_shared:,} regions")
print(f"  With CLOSED human orthologs: {mouse_specific:,} regions")
print(f"  Total classified:            {mouse_shared + mouse_specific:,} regions")

#Visualization: Open vs Closed orthologs
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

#Human perspective
axes[0].pie(
    [human_shared, human_specific],
    labels=[f"Open in Mouse\n({human_shared:,})", f"Closed in Mouse\n({human_specific:,})"],
    autopct="%1.1f%%",  #Show percentages on pie chart
    colors=["#99ff99", "#ff9999"],
    startangle=90
)
axes[0].set_title("Human OCRs: Ortholog Status in Mouse", fontsize=12, fontweight="bold")

#Mouse perspective
axes[1].pie(
    [mouse_shared, mouse_specific],
    labels=[f"Open in Human\n({mouse_shared:,})", f"Closed in Human\n({mouse_specific:,})"],
    autopct="%1.1f%%",
    colors=["#99ff99", "#ff9999"],
    startangle=90
)
axes[1].set_title("Mouse OCRs: Ortholog Status in Human", fontsize=12, fontweight="bold")

plt.tight_layout()
plt.savefig("ortholog_status.png", dpi=300, bbox_inches="tight")
print("Saved: ortholog_status.png")


#PROMOTER AND ENHANCER CLASSIFICATION
print("\nPromoter and Enhancer Classification")

#Read all classification files
human_prom = len(pd.read_csv(results_dir / "human_all_promoters.bed", sep="\t", header=None))
human_enh = len(pd.read_csv(results_dir / "human_all_enhancers.bed", sep="\t", header=None))
human_total = human_prom + human_enh

mouse_prom = len(pd.read_csv(results_dir / "mouse_all_promoters.bed", sep="\t", header=None))
mouse_enh = len(pd.read_csv(results_dir / "mouse_all_enhancers.bed", sep="\t", header=None))
mouse_total = mouse_prom + mouse_enh

shared_prom = len(pd.read_csv(results_dir / "shared_promoters.bed", sep="\t", header=None))
shared_enh = len(pd.read_csv(results_dir / "shared_enhancers.bed", sep="\t", header=None))
shared_total = shared_prom + shared_enh

human_spec_prom = len(pd.read_csv(results_dir / "human_specific_promoters.bed", sep="\t", header=None))
human_spec_enh = len(pd.read_csv(results_dir / "human_specific_enhancers.bed", sep="\t", header=None))
human_spec_total = human_spec_prom + human_spec_enh

mouse_spec_prom = len(pd.read_csv(results_dir / "mouse_specific_promoters.bed", sep="\t", header=None))
mouse_spec_enh = len(pd.read_csv(results_dir / "mouse_specific_enhancers.bed", sep="\t", header=None))
mouse_spec_total = mouse_spec_prom + mouse_spec_enh

print("\nClassification Method:")
print("  - Promoters: OCRs within 2000 bp of nearest TSS")
print("  - Enhancers: OCRs beyond 2000 bp from nearest TSS")
print("  - Method: Distance-based using bedtools closest")
print("\nLimitations:")
print("  - Distance threshold is arbitrary (2kb is commonly used but not universal)")
print("  - Does not account for gene density or chromatin state")
print("  - Some true promoters may extend beyond 2kb")
print("  - Some enhancers may be closer than 2kb to TSS")

print("\nAll OCRs by species:")
print(f"  Human: {human_prom:,} promoters ({100*human_prom/human_total:.1f}%), "
      f"{human_enh:,} enhancers ({100*human_enh/human_total:.1f}%)")
print(f"  Mouse: {mouse_prom:,} promoters ({100*mouse_prom/mouse_total:.1f}%), "
      f"{mouse_enh:,} enhancers ({100*mouse_enh/mouse_total:.1f}%)")

print("\nShared (conserved) OCRs:")
print(f"  {shared_prom:,} promoters ({100*shared_prom/shared_total:.1f}%), "
      f"{shared_enh:,} enhancers ({100*shared_enh/shared_total:.1f}%)")

print("\nSpecies-specific OCRs:")
print(f"  Human: {human_spec_prom:,} promoters ({100*human_spec_prom/human_spec_total:.1f}%), "
      f"{human_spec_enh:,} enhancers ({100*human_spec_enh/human_spec_total:.1f}%)")
print(f"  Mouse: {mouse_spec_prom:,} promoters ({100*mouse_spec_prom/mouse_spec_total:.1f}%), "
      f"{mouse_spec_enh:,} enhancers ({100*mouse_spec_enh/mouse_spec_total:.1f}%)")

#Visualization: Promoter vs Enhancer ratios
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

#All OCRs
axes[0].pie(
    [human_prom, human_enh],
    labels=[f"Promoters\n({human_prom:,})", f"Enhancers\n({human_enh:,})"],
    autopct="%1.1f%%",  #Show percentages on pie chart
    colors=["#ff9999", "#66b3ff"],
    startangle=90
)
axes[0].set_title(f"Human All OCRs (Total: {human_total:,})", fontsize=12, fontweight="bold")

axes[1].pie(
    [mouse_prom, mouse_enh],
    labels=[f"Promoters\n({mouse_prom:,})", f"Enhancers\n({mouse_enh:,})"],
    autopct="%1.1f%%",
    colors=["#ff9999", "#66b3ff"],
    startangle=90
)
axes[1].set_title(f"Mouse All OCRs (Total: {mouse_total:,})", fontsize=12, fontweight="bold")

plt.tight_layout()
plt.savefig("promoter_enhancer_classification.png", dpi=300, bbox_inches="tight")
print("Saved: promoter_enhancer_classification.png")


#CONSERVATION ANALYSIS
print("\nConservation Analysis")

print(f"\nShared (conserved) OCRs: {shared_total:,} total")
print(f"  Promoters: {shared_prom:,} ({100*shared_prom/shared_total:.2f}%)")
print(f"  Enhancers: {shared_enh:,} ({100*shared_enh/shared_total:.2f}%)")
print(f"  Enhancer/Promoter Ratio: {shared_enh/shared_prom:.2f}")

print(f"\nHuman-specific OCRs: {human_spec_total:,} total")
print(f"  Promoters: {human_spec_prom:,} ({100*human_spec_prom/human_spec_total:.2f}%)")
print(f"  Enhancers: {human_spec_enh:,} ({100*human_spec_enh/human_spec_total:.2f}%)")
print(f"  Enhancer/Promoter Ratio: {human_spec_enh/human_spec_prom:.2f}")

print(f"\nMouse-specific OCRs: {mouse_spec_total:,} total")
print(f"  Promoters: {mouse_spec_prom:,} ({100*mouse_spec_prom/mouse_spec_total:.2f}%)")
print(f"  Enhancers: {mouse_spec_enh:,} ({100*mouse_spec_enh/mouse_spec_total:.2f}%)")
print(f"  Enhancer/Promoter Ratio: {mouse_spec_enh/mouse_spec_prom:.2f}")

#Visualization: Conservation by element type
fig, ax = plt.subplots(figsize=(10, 6))

categories = ['Shared\n(Conserved)', 'Human\nSpecific', 'Mouse\nSpecific']
promoters = [shared_prom, human_spec_prom, mouse_spec_prom]
enhancers = [shared_enh, human_spec_enh, mouse_spec_enh]

x = np.arange(len(categories))
width = 0.6

p1 = ax.bar(x, promoters, width, label='Promoters', color='#ff9999')
p2 = ax.bar(x, enhancers, width, bottom=promoters, label='Enhancers', color='#66b3ff')

ax.set_ylabel('Number of Regions', fontsize=12)
ax.set_title('Conservation Status: Promoters vs Enhancers', fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.legend()

#Add total labels on bars
for i, (p, e) in enumerate(zip(promoters, enhancers)):
    total = p + e
    ax.text(i, total + 50, f'{total:,}', ha='center', fontweight='bold')

plt.tight_layout()
plt.savefig("conservation_comparison.png", dpi=300, bbox_inches="tight")
print("Saved: conservation_comparison.png")

#Additional visualization: Percentage comparison
fig, ax = plt.subplots(figsize=(10, 6))

prom_pcts = [
    100 * shared_prom / shared_total,
    100 * human_spec_prom / human_spec_total,
    100 * mouse_spec_prom / mouse_spec_total
]
enh_pcts = [
    100 * shared_enh / shared_total,
    100 * human_spec_enh / human_spec_total,
    100 * mouse_spec_enh / mouse_spec_total
]

x = np.arange(len(categories))
width = 0.35

bars1 = ax.bar(x - width/2, prom_pcts, width, label='Promoters', color='#ff9999')
bars2 = ax.bar(x + width/2, enh_pcts, width, label='Enhancers', color='#66b3ff')

ax.set_ylabel('Percentage (%)', fontsize=12)
ax.set_title('Promoter vs Enhancer Distribution by Conservation Status', 
             fontsize=14, fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.legend()
ax.set_ylim(0, 100)

#Add percentage labels on bars
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{height:.1f}%', ha='center', va='bottom', fontsize=10)

plt.tight_layout()
plt.savefig("percentage_comparison.png", dpi=300, bbox_inches="tight")
print("Saved: percentage_comparison.png")


#SUMMARY TABLE
print("\nSUMMARY TABLE")

#Create structured data for summary table
summary_data = {
    "Category": [
        "Human All",
        "Mouse All",
        "Shared (Conserved)",
        "Human Specific",
        "Mouse Specific",
    ],
    "Promoters": [
        human_prom,
        mouse_prom,
        shared_prom,
        human_spec_prom,
        mouse_spec_prom,
    ],
    "Enhancers": [
        human_enh,
        mouse_enh,
        shared_enh,
        human_spec_enh,
        mouse_spec_enh,
    ],
    "Total": [
        human_total,
        mouse_total,
        shared_total,
        human_spec_total,
        mouse_spec_total,
    ],
    "Promoter %": [
        100 * human_prom / human_total,
        100 * mouse_prom / mouse_total,
        100 * shared_prom / shared_total,
        100 * human_spec_prom / human_spec_total,
        100 * mouse_spec_prom / mouse_spec_total,
    ],
    "Enhancer %": [
        100 * human_enh / human_total,
        100 * mouse_enh / mouse_total,
        100 * shared_enh / shared_total,
        100 * human_spec_enh / human_spec_total,
        100 * mouse_spec_enh / mouse_spec_total,
    ]
}

summary_df = pd.DataFrame(summary_data)
summary_df = summary_df.round(2)
print(summary_df.to_string(index=False))

#Save to CSV
summary_df.to_csv("summary_table.csv", index=False)
print("Saved: summary_table.csv")


#SANITY CHECKS
print("\nSANITY CHECKS")

checks_passed = 0
checks_total = 0

def check(name, condition, details=""):
    global checks_passed, checks_total
    checks_total += 1
    symbol = "✓" if condition else "✗"
    status = "PASS" if condition else "FAIL"
    print(f"{symbol} [{status}] {name}")
    if details and not condition:
        print(f"         {details}")
    if condition:
        checks_passed += 1

#Check 1: Promoter and enhancer counts should sum to total
check("Human promoters + enhancers = total",
      human_prom + human_enh == human_total)

check("Mouse promoters + enhancers = total",
      mouse_prom + mouse_enh == mouse_total)

check("Shared promoters + enhancers = total",
      shared_prom + shared_enh == shared_total)

#Check 2: Promoter to enhancer ratio should be reasonable (10-50% promoters typical)
human_prom_pct = 100 * human_prom / human_total
check("Human promoter % reasonable (10-50%)",
      10 <= human_prom_pct <= 50,
      f"Got {human_prom_pct:.1f}%")

mouse_prom_pct = 100 * mouse_prom / mouse_total
check("Mouse promoter % reasonable (10-50%)",
      10 <= mouse_prom_pct <= 50,
      f"Got {mouse_prom_pct:.1f}%")

#Check 3: There should be some conservation
check("Shared regions exist", shared_total > 0)

#Check 4: There should be species-specific regions
check("Human-specific regions exist", human_spec_total > 0)
check("Mouse-specific regions exist", mouse_spec_total > 0)

#Check 5: Conservation should not be 100%
check("Not all regions are conserved",
      shared_total < min(human_total, mouse_total))

print(f"\nChecks passed: {checks_passed}/{checks_total}")


#FINISHED
print("\nAnalysis Complete")
print("\nGenerated files:")
print("  1. ortholog_status.png")
print("  2. promoter_enhancer_classification.png")
print("  3. conservation_comparison.png")
print("  4. percentage_comparison.png")
print("  5. summary_table.csv")