#!/usr/bin/env python3

from pathlib import Path
import argparse
import subprocess
import sys

DEFAULT_BEDS = [
    "human_promoters.bed",
    "human_enhancers.bed",
    "human_mapped_from_mouse_promoters.bed",
    "human_mapped_from_mouse_enhancers.bed",
]

def main():
    parser = argparse.ArgumentParser(
        description="Batch-run GREAT analysis on BED files using rGREAT online mode."
    )
    parser.add_argument(
        "--bed-dir",
        default="results/classification_results/raw_results",
        help="Directory containing BED files."
    )
    parser.add_argument(
        "--outdir",
        default="enrichment_analysis/great_results",
        help="Directory to store GREAT results."
    )
    parser.add_argument(
        "--species",
        default="hg38",
        help="Genome assembly/species for GREAT (default: hg38)."
    )
    parser.add_argument(
        "--beds",
        nargs="+",
        default=DEFAULT_BEDS,
        help="BED filenames to process."
    )
    parser.add_argument(
        "--rscript-bin",
        default="Rscript",
        help="Path to Rscript executable."
    )
    parser.add_argument(
        "--great-r-script",
        default="enrichment_analysis/run_great_online.R",
        help="Path to R script that runs online GREAT via rGREAT."
    )

    args = parser.parse_args()

    bed_dir = Path(args.bed_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for bed_name in args.beds:
        bed_path = bed_dir / bed_name
        if not bed_path.exists():
            print(f"[WARN] Missing BED file, skipping: {bed_path}", file=sys.stderr)
            continue

        prefix = bed_path.stem
        sample_outdir = outdir / prefix
        sample_outdir.mkdir(parents=True, exist_ok=True)

        cmd = [
            args.rscript_bin,
            args.great_r_script,
            "--bed", str(bed_path),
            "--species", args.species,
            "--outdir", str(sample_outdir),
            "--prefix", prefix,
        ]

        print(f"[INFO] Running GREAT for {bed_path}", file=sys.stderr)
        subprocess.run(cmd, check=True)

    print("[INFO] GREAT batch complete.", file=sys.stderr)

if __name__ == "__main__":
    main()