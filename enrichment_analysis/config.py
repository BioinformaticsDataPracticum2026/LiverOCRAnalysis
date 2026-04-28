from pathlib import Path

# folder where the input BED files (from classification step) are stored
INPUT_BED_DIR = Path("results/classification_results/raw_results")

GREAT_OUTPUT_DIR = Path("results/enrichment/great")

# mapping of dataset names to their BED file and genome
BED_FILES = {
    # human OCRs that are shared between species
    "human_shared": {"file": "human_shared.bed", "genome": "hg38"},

    # human OCRs that are species-specific
    "human_specific": {"file": "human_specific.bed", "genome": "hg38"},

    # mouse OCRs that are shared between species
    "mouse_shared": {"file": "mouse_shared.bed", "genome": "mm10"},

    # mouse OCRs that are species-specific
    "mouse_specific": {"file": "mouse_specific.bed", "genome": "mm10"},

    # shared OCRs classified as promoters
    "shared_promoters": {"file": "shared_promoters.bed", "genome": "mm10"},

    # shared OCRs classified as enhancers
    "shared_enhancers": {"file": "shared_enhancers.bed", "genome": "mm10"},
}