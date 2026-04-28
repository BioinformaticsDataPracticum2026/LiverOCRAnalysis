import subprocess
from pathlib import Path
from enrichment_analysis.config import INPUT_BED_DIR, GREAT_OUTPUT_DIR, BED_FILES


# path to the R script that runs GREAT
R_SCRIPT = "enrichment_analysis/great_online.R"


def run_one(label, bed_path, genome):
    # create output folder for this specific dataset
    outdir = GREAT_OUTPUT_DIR / label
    outdir.mkdir(parents=True, exist_ok=True)

    # command to run the R script with inputs
    cmd = [
        "Rscript",
        R_SCRIPT,
        str(bed_path),
        genome,
        str(outdir)
    ]

    # print which dataset is currently running
    print(f"[RUN] {label}")

    # run the command in the terminal
    subprocess.run(cmd, check=True)


def main():
    # loop through all BED files defined in config
    for label, info in BED_FILES.items():
        bed_path = INPUT_BED_DIR / info["file"]  # path to the BED file
        genome = info["genome"]  # genome version (e.g., hg38)

        # run GREAT for this dataset
        run_one(label, bed_path, genome)


if __name__ == "__main__":
    main()