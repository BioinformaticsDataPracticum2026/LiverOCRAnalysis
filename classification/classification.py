#!/usr/bin/env python3

"""
Pipeline poriton to classify OCRs as promoters/enhancers and compare regulatory states across species.

Usage:
    ./classify_ocr.py --config /path/to/config.yaml
    
SLURM Script:
    #!/bin/bash
    #SBATCH -J classification
    #SBATCH -p RM-shared
    #SBATCH -N 1
    #SBATCH -n 2
    #SBATCH -t 12:00:00
    #SBATCH --mem-per-cpu=2000M 
    #SBATCH -o logs/classification_%j.out 
    #SBATCH -e logs/classification_%j.err
    #SBATCH --mail-type=END,FAIL
    #SBATCH --mail-user=alfredl@andrew.cmu.edu
    module load anaconda3
    conda activate my_bio_env
    ./classify_ocr.py --config /path/to/config.yaml
"""

from pathlib import Path
from pybedtools import BedTool #Available in bioconda environment
from main import run_classification
from typing import Dict, Tuple, Optional
import logging
import argparse

#Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)
logger = logging.getLogger(__name__)

#Helper Functions
def classifyOcrPromotersEnhancers(
    ocr_bed_path: str,
    tss_bed_path: str,
    output_prefix: str,
    promoter_distance: int = 2000  #Might be worth switching to 5k?
) -> Dict[str, str]:
    """
    Classify OCRs as promoters or enhancers based on distance to nearest TSS.
    Promoters defined as OCRs within promoter_distance bp of a TSS.
    Enhancers are OCRs beyond that distance.
    
    Inputs:
        ocr_bed_path : str -> Path to the OCR BED file (peaks).
        tss_bed_path : str -> Path to the TSS BED file (1bp per TSS).
        output_prefix : str -> Prefix for output files. Outputs will be saved as:
            - {output_prefix}_promoters.bed
            - {output_prefix}_enhancers.bed
            - {output_prefix}_nearestTSS.bed
        promoter_distance : int, optional -> Distance threshold for promoter classification (default: 2000 bp).
            Consider 5000 bp for more lenient classification.
    
    Outputs:
        Dict[str, str] -> Dictionary with keys 'promoters', 'enhancers', 'nearest' containing paths to the output BED files.
    
    Errors:
        FileNotFoundError -> If input BED files do not exist.
        ValueError -> If BED files are empty or invalid.
    """
    
    ocr_bed_path = Path(ocr_bed_path)
    tss_bed_path = Path(tss_bed_path)

    #Validate input files
    if not ocr_bed_path.exists():
        raise FileNotFoundError(f"OCR BED file not found: {ocr_bed_path}")
    if not tss_bed_path.exists():
        raise FileNotFoundError(f"TSS BED file not found: {tss_bed_path}")
    
    logging.info(f"Processing OCRs: {ocr_bed_path}")
    
    try:
        #Load BED files
        ocr = BedTool(ocr_bed_path)
        tss = BedTool(tss_bed_path)

        #Validate non-empty
        if len(ocr) == 0:
            raise ValueError("OCR BED file is empty")
        if len(tss) == 0:
            raise ValueError("TSS BED file is empty")
        
        logger.info(f"Loaded {len(ocr)} OCRs and {len(tss)} TSSs")

        #Annotate nearest TSS + distance
        ocr_with_genes = ocr.closest(tss, d=True)
        nearest_tss_path = f"{output_prefix}_nearestTSS.bed"
        ocr_with_genes.saveas(nearest_tss_path)
        logger.info(f"Saved annotated OCRs with nearest TSS: {nearest_tss_path}")
    
        #Classify OCRs into promoters by distance
        promoters_path = f"{output_prefix}_promoters.bed"
        ocr_promoters = ocr_with_genes.filter(lambda x: int(x[-1]) <= promoter_distance).saveas(promoters_path)
        logger.info(
            f"Classified {len(ocr_promoters)} promoters (<=  {promoter_distance} bp): "
            f"{promoters_path}"
        )

        #Classify OCRs into enhancers by distance
        enhancers_path = f"{output_prefix}_enhancers.bed"
        ocr_enhancers = ocr_with_genes.filter(lambda x: int(x[-1]) > promoter_distance).saveas(f"{output_prefix}_enhancers.bed")
        logger.info(
            f"Classified {len(ocr_enhancers)} enhancers (> {promoter_distance} bp): "
            f"{enhancers_path}"
        )

        logging.info(f"Finished processing {ocr_bed_path}")


        return {
            "promoters": promoters_path,
            "enhancers": enhancers_path,
            "nearest": nearest_tss_path
        }
    
    except Exception as e:
        logger.error(f"Error processing {ocr_bed_path}: {e}")
        raise

#classifyConservedRegions() --> Assigns promoter/enhancer labels to mapped (conserved) regions using TSS distance.
def classifyConservedRegions(
        conservedBedPath: str,
        tssBedPath: str,
        outputPrefix: str,
        promoterDistance: int = 2000
        ):

    logging.info(f"Processing conserved regions: {conservedBedPath}")

    conserved = BedTool(conservedBedPath)
    tss = BedTool(tssBedPath)

    # Annotate distance
    annotated = conserved.closest(tss, d=True)
    annotated.saveas(f"{outputPrefix}_TSS.bed")

    # Split promoter/enhancer
    promoters = annotated.filter(lambda x: int(x[-1]) <= promoterDistance).saveas(
        f"{outputPrefix}_promoters.bed"
    )
    enhancers = annotated.filter(lambda x: int(x[-1]) >  promoterDistance).saveas(
        f"{outputPrefix}_enhancers.bed"
    )

    return promoters, enhancers

#findSharedElements() --> Identify conserved regions that overlap native regulatory elements in the target species.
def findSharedElements(
        mappedFilePath: str,
        nativeFilePath: str,
        outputFilePath: str
        ):

    logging.info(f"Finding shared elements")

    mapped = BedTool(mappedFilePath)
    native = BedTool(nativeFilePath)

    shared = mapped.intersect(native, u=True)
    shared.saveas(outputFilePath)

    return outputFilePath

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Classify OCRs into promoters/enhancers")
    parser.add_argument("--config", required=True)

    args = parser.parse_args()

    run_classification(Path(args.config))
