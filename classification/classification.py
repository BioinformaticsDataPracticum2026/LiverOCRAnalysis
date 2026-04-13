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

def classifyConservedRegions(
        conserved_bed_path: str,
        tss_bed_path: str,
        output_prefix: str,
        promoter_distance: int = 2000
) -> Tuple[str, str]:
    """
    Classify conserved/mapped regions as promoters or enhancers using TSS distance.
    
    Inputs:
        conserved_bed_path : str -> Path to the conserved/mapped OCR BED file.
        tss_bed_path : str -> Path to the TSS BED file (1bp per TSS).
        output_prefix : str -> Prefix for output files. Outputs will be saved as:
            - {output_prefix}_promoters.bed
            - {output_prefix}_enhancers.bed
            - {output_prefix}_TSS.bed (annotated)
        promoter_distance : int, optional -> Distance threshold for promoter classification (default: 2000 bp).
    
    Outputs:
        Tuple[str, str] -> Paths to promoters and enhancers BED files.
    
    Errors:
        FileNotFoundError -> If input BED files do not exist.
        ValueError -> If BED files are empty or invalid.
    """
    conserved_bed_path = Path(conserved_bed_path)
    tss_bed_path = Path(tss_bed_path)

    if not conserved_bed_path.exists():
        raise FileNotFoundError(f"Conserved BED file not found: {conserved_bed_path}")
    if not tss_bed_path.exists():
        raise FileNotFoundError(f"TSS BED file not found: {tss_bed_path}")

    logging.info(f"Processing conserved regions: {conserved_bed_path}")

    try:
        conserved = BedTool(conserved_bed_path)
        tss = BedTool(tss_bed_path)

        if len(conserved) == 0:
            raise ValueError("Conserved BED file is empty")
        if len(tss) == 0:
            raise ValueError("TSS BED file is empty")
        
        logger.info(f"Loaded {len(conserved)} conserved regions and {len(tss)} TSSs")

        #Annotate distance to nearest TSS
        tss_path = f"{output_prefix}_TSS.bed"
        annotated = conserved.closest(tss, d=True)
        annotated.saveas(tss_path)
        logger.info(f"Saved annotated conserved regions: {tss_path}")

        #Split into promoters
        promoters_path = f"{output_prefix}_promoters.bed"
        promoters = annotated.filter(lambda x: int(x[-1]) <= promoter_distance).saveas(promoters_path)
        logger.info(
            f"Classified {len(promoters)} conserved promoters "
            f"(<= {promoter_distance} bp): {promoters_path}"
        )

        #Split into enhancers
        enhancers_path = f"{output_prefix}_enhancers.bed"
        enhancers = annotated.filter(lambda x: int(x[-1]) >  promoter_distance).saveas(enhancers_path)
        logger.info(
            f"Classified {len(enhancers)} conserved enhancers "
            f"(> {promoter_distance} bp): {enhancers_path}"
        )

        return promoters_path, enhancers_path
    
    except Exception as e:
        logger.error(f"Error processing {conserved_bed_path}: {e}")
        raise


def findSharedElements(
        mapped_file_path: str,
        native_file_path: str,
        output_file_path: str
) -> str:
    """
    Identify conserved regions that overlap native regulatory elements in target species.
        - Finds mapped OCRs from one species that intersect with native
        OCRs detected in another species, identifying true conserved regulatory elements.
    
    Inputs:
        mapped_file_path : str -> Path to mapped/conserved OCR BED file (from source species, mapped to target).
        native_file_path : str -> Path to native OCR BED file (detected in target species).
        output_file_path : str -> Path for output BED file containing overlapping regions.
    
    Ouputs:
        str -> Path to the output BED file with shared elements.
    
    Errors:
        FileNotFoundError -> If input BED files do not exist.
        ValueError -> If BED files are empty or invalid.
    """
    mapped_file_path = Path(mapped_file_path)
    native_file_path = Path(native_file_path)

    if not mapped_file_path.exists():
        raise FileNotFoundError(f"Mapped BED file not found: {mapped_file_path}")
    if not native_file_path.exists():
        raise FileNotFoundError(f"Native BED file not found: {native_file_path}")

    logging.info(f"Finding shared elements between {mapped_file_path} and {native_file_path}")

    try:
        mapped = BedTool(mapped_file_path)
        native = BedTool(native_file_path)

        if len(mapped) == 0:
            logger.warning("Mapped BED file is empty")
        if len(native) == 0:
            logger.warning("Native BED file is empty")        

        logger.info(f"Loaded {len(mapped)} mapped regions and {len(native)} native regions")

        #Find intersection regions
        shared = mapped.intersect(native, u=True) #u=True -> Keeps only regions that intersect
        shared.saveas(output_file_path)

        logger.info(f"Found {len(shared)} shared elements: {output_file_path}")

        return output_file_path

    except Exception as e:
        logger.error(f"Error finding shared elements: {e}")
        raise

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Classify OCRs into promoters/enhancers")
    parser.add_argument("--config", required=True)

    args = parser.parse_args()

    run_classification(Path(args.config))
