import subprocess
from pathlib import Path
from enrichment_analysis.config import INPUT_BED_DIR, GREAT_OUTPUT_DIR, BED_FILES


R_SCRIPT = "enrichment_analysis/great_online.R"


def run_one(label, bed_path, genome):
    outdir = GREAT_OUTPUT_DIR / label
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        "Rscript",
        R_SCRIPT,
        str(bed_path),
        genome,
        str(outdir)
    ]

    print(f"[RUN] {label}")
    subprocess.run(cmd, check=True)


def main():
    for label, info in BED_FILES.items():
        bed_path = INPUT_BED_DIR / info["file"]
        genome = info["genome"]

        run_one(label, bed_path, genome)


if __name__ == "__main__":
    main()