#!/usr/bin/env python3

"""
Pipeline poriton to classify OCRs as promoters/enhancers and compare regulatory states across species.

Usage:
    ./classify_ocr.py --config /path/to/config.yaml
    
Sample SLURM Script:
    #!/bin/bash
    #SBATCH -J ocr_classification
    #SBATCH -p RM-shared
    #SBATCH -N 1
    #SBATCH -n 2
    #SBATCH -t 12:00:00
    #SBATCH --mem-per-cpu=2000M 
    #SBATCH -o logs/classification_%j.out 
    #SBATCH -e logs/classification_%j.err
    #SBATCH --mail-type=END,FAIL
    #SBATCH --mail-user=alfredl@andrew.cmu.edu

    #Load anaconda and activate environment
    module load anaconda3
    conda activate my_bio_env

    #Run the classification pipeline
    python3 classification.py --config config.processed.yaml
"""

from pathlib import Path
from pybedtools import BedTool #Available in bioconda environment
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
    promoter_distance: int = 2000,
    save_nearest: bool = False
) -> Dict[str, str]:
    """
    Classify OCRs as promoters or enhancers based on distance to nearest TSS.
    Promoters are OCRs within promoter_distance of a TSS.
    Enhancers are OCRs beyond that distance.
    
    Inputs:
        ocr_bed_path : str -> Path to the OCR BED file (peaks).
        tss_bed_path : str -> Path to the TSS BED file (1bp per TSS).
        output_prefix : str -> Prefix for output files. Outputs will be saved as:
            - {output_prefix}_promoters.bed
            - {output_prefix}_enhancers.bed
            - {output_prefix}_nearestTSS.bed (only if save_nearest=True)
        promoter_distance : int, optional -> Distance threshold for promoter classification (default: 2000 bp).
        save_nearest : bool, optional -> Whether to save intermediate nearestTSS file (default: False).
    
    Outputs:
        Dict[str, str] -> Dictionary with keys 'promoters' and 'enhancers' containing paths to the output BED files.
    
    Errors:
        FileNotFoundError -> If input BED files do not exist.
        ValueError -> If BED files are empty or invalid.
    """
    
    #Assign file paths to variables
    ocr_bed_path = Path(ocr_bed_path)
    tss_bed_path = Path(tss_bed_path)

    #Validate input files exist
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

        # Only save nearestTSS file if requested (for debugging)
        if save_nearest:
            nearest_tss_path = f"{output_prefix}_nearestTSS.bed"
            ocr_with_genes.saveas(nearest_tss_path)
            logger.info(f"Saved annotated OCRs with nearest TSS: {nearest_tss_path}")

    
        #Classify OCRs into promoters by distance
        promoters_path = f"{output_prefix}_promoters.bed"
        ocr_promoters = ocr_with_genes.filter(lambda x: int(x[-1]) <= promoter_distance).saveas(promoters_path)
        logger.info(
            f"Classified {len(ocr_promoters)} promoters (<=  {promoter_distance} bp): {promoters_path}"
        )

        #Classify OCRs into enhancers by distance
        enhancers_path = f"{output_prefix}_enhancers.bed"
        ocr_enhancers = ocr_with_genes.filter(lambda x: int(x[-1]) > promoter_distance).saveas(f"{output_prefix}_enhancers.bed")
        logger.info(
            f"Classified {len(ocr_enhancers)} enhancers (> {promoter_distance} bp): {enhancers_path}"
        )

        logging.info(f"Finished processing {ocr_bed_path}")

        return {
            "promoters": promoters_path,
            "enhancers": enhancers_path,
        }
    
    #Error if promoter/enhancer classification fails
    except Exception as e:
        logger.error(f"Error processing {ocr_bed_path}: {e}")
        raise

def identifyOrthologStatus(
    source_ocr_bed_path: str,
    source_to_target_mapped_path: str,
    target_native_ocr_path: str,
    output_prefix: str
) -> Dict[str, str]:
    """
    Classify source species OCRs by ortholog accessibility status in target species.
    
    Identifies:
    - Shared: Source OCRs whose orthologs are OPEN in target species
    - Specific: Source OCRs whose orthologs are CLOSED in target species
    
    Inputs:
        source_ocr_bed_path : str -> Path to source species OCR BED file.
        source_to_target_mapped_path : str -> Path to HALPER mapped regions.
        target_native_ocr_path : str -> Path to target species native OCR BED file.
        output_prefix : str -> Prefix for output files. Outputs will be saved as:
            - {output_prefix}_shared.bed (open in both species)
            - {output_prefix}_specific.bed (open in source, closed in target)
    
    Outputs:
        Dict[str, str] -> Dictionary with keys 'shared' and 'specific' containing output paths.
    
    Errors:
        FileNotFoundError -> If input BED files do not exist.
        ValueError -> If BED files are empty or invalid.
    """
    
    #Assign file paths to variables
    source_ocr_bed_path = Path(source_ocr_bed_path)
    source_to_target_mapped_path = Path(source_to_target_mapped_path)
    target_native_ocr_path = Path(target_native_ocr_path)
    
    #Validate input files exist
    if not source_ocr_bed_path.exists():
        raise FileNotFoundError(f"Source OCR BED file not found: {source_ocr_bed_path}")
    if not source_to_target_mapped_path.exists():
        raise FileNotFoundError(f"Mapped BED file not found: {source_to_target_mapped_path}")
    if not target_native_ocr_path.exists():
        raise FileNotFoundError(f"Target OCR BED file not found: {target_native_ocr_path}")
    
    logger.info(f"Identifying ortholog status for {source_ocr_bed_path}")
    
    try:
        #Load BED files
        source_ocr = BedTool(source_ocr_bed_path)
        mapped = BedTool(source_to_target_mapped_path)
        target_native = BedTool(target_native_ocr_path)
        
        #Validate non-empty
        if len(source_ocr) == 0:
            raise ValueError("Source OCR BED file is empty")
        if len(mapped) == 0:
            logger.warning("Mapped BED file is empty")
        if len(target_native) == 0:
            logger.warning("Target OCR BED file is empty")
        
        logger.info(
            f"Loaded {len(source_ocr)} source OCRs, "
            f"{len(mapped)} mapped regions, "
            f"{len(target_native)} target OCRs"
        )
        
        #Find mapped regions that overlap target OCRs (open in both = shared)
        mapped_open = mapped.intersect(target_native, u=True)
        
        #Find mapped regions that DON'T overlap target OCRs (closed in target = specific)
        mapped_closed = mapped.intersect(target_native, v=True)
        
        logger.info(
            f"Mapped regions: {len(mapped_open)} open in target, "
            f"{len(mapped_closed)} closed in target"
        )
        
        #Trace back to original source coordinates
        source_shared = source_ocr.intersect(mapped_open, u=True)
        source_specific = source_ocr.intersect(mapped_closed, u=True)
        
        #Save outputs
        shared_path = f"{output_prefix}_shared.bed"
        specific_path = f"{output_prefix}_specific.bed"
        
        source_shared.saveas(shared_path)
        source_specific.saveas(specific_path)
        
        logger.info(f"Saved {len(source_shared)} shared OCRs: {shared_path}")
        logger.info(f"Saved {len(source_specific)} species-specific OCRs: {specific_path}")
        
        return {
            "shared": shared_path,
            "specific": specific_path
        }
    
    except Exception as e:
        logger.error(f"Error identifying ortholog status: {e}")
        raise


#Main wrapper function
def run_classification(config_path: Path) -> None:
    """
    Main pipeline orchestrator. Loads config and runs classification steps.
    
    Produces 10 essential BED files:
    - {species}_all_promoters.bed
    - {species}_all_enhancers.bed
    - {species}_specific_promoters.bed
    - {species}_specific_enhancers.bed
    - shared_promoters.bed
    - shared_enhancers.bed
    
    Plus summary_stats.txt answering Goal 4a.
    
    Expected YAML config structure:
    ```
    species_1: "human"
    species_2: "mouse"
    
    species_1_peak_file_cleaned: path/to/species1_ocrs.bed
    species_1_tss_file_cleaned: path/to/species1_tss.bed
    species_2_peak_file_cleaned: path/to/species2_ocrs.bed
    species_2_tss_file_cleaned: path/to/species2_tss.bed
    
    species_1_to_species_2_cleaned: path/to/species1_to_species2_mapped.bed
    species_2_to_species_1_cleaned: path/to/species2_to_species1_mapped.bed
    
    output_dir: "classification_results"
    
    parameters:
      promoter_distance: 2000
    ```
    
    Inputs:
        config_path : Path -> Path to YAML configuration file.
    """
    #Import yaml -> Config is a yaml file
    try:
        import yaml
    except ImportError:
        logger.error("PyYAML not found. Install with: pip install pyyaml")
        raise
    
    #Assign config file path to variable
    config_path = Path(config_path)

    #Validate config file exists
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    logger.info(f"Loading configuration from {config_path}")
    
    #Read in config file
    try:
        with open(config_path) as f:
            config = yaml.safe_load(f)
    except Exception as e:
        logger.error(f"Error loading config: {e}")
        raise
    
    #Extract parameters
    params = config.get("parameters", {})
    promoter_distance = params.get("promoter_distance", 2000)
    output_dir = Path(config.get("output_dir", "classification_results"))
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Promoter distance threshold: {promoter_distance} bp")

    #Build species config
    #Get species names
    species_1_name = config.get("species_1")
    species_2_name = config.get("species_2")

    logger.info(f"Species 1: {species_1_name}")
    logger.info(f"Species 2: {species_2_name}")

    #Step 1: Classify all OCRS by type
    logger.info("Step 1: Classifying all OCRs by type")

    s1_results = classifyOcrPromotersEnhancers(
        ocr_bed_path=config.get("species_1_peak_file_cleaned"),
        tss_bed_path=config.get("species_1_tss_file_cleaned"),
        output_prefix=str(output_dir / f"{species_1_name}_all"),
        promoter_distance=promoter_distance,
        save_nearest=False  # Don't save intermediate files
    )
    
    s2_results = classifyOcrPromotersEnhancers(
        ocr_bed_path=config.get("species_2_peak_file_cleaned"),
        tss_bed_path=config.get("species_2_tss_file_cleaned"),
        output_prefix=str(output_dir / f"{species_2_name}_all"),
        promoter_distance=promoter_distance,
        save_nearest=False
    )

    #Organizing results for logging
    classification_results = {
        species_1_name: {
            "all_promoters": s1_results["promoters"],
            "all_enhancers": s1_results["enhancers"]
        },
        species_2_name: {
            "all_promoters": s2_results["promoters"],
            "all_enhancers": s2_results["enhancers"]
        }
    }


    #Step 2: Identify shared and species-specific OCRs
    logger.info("Step 2: Identifying shared and species-specific OCRs")
    
    #Species 1 perspective: which are shared vs specific?
    s1_ortholog_status = identifyOrthologStatus(
        source_ocr_bed_path=config.get("species_1_peak_file_cleaned"),
        source_to_target_mapped_path=config.get("species_1_to_species_2_cleaned"),
        target_native_ocr_path=config.get("species_2_peak_file_cleaned"),
        output_prefix=str(output_dir / f"{species_1_name}")
    )
    
    #Species 2 perspective: which are shared vs specific?
    s2_ortholog_status = identifyOrthologStatus(
        source_ocr_bed_path=config.get("species_2_peak_file_cleaned"),
        source_to_target_mapped_path=config.get("species_2_to_species_1_cleaned"),
        target_native_ocr_path=config.get("species_1_peak_file_cleaned"),
        output_prefix=str(output_dir / f"{species_2_name}")
    )

    #Organizing results for logging
    ortholog_results = {
        f"{species_1_name}_ortholog_status": {
            "shared": s1_ortholog_status["shared"],
            "specific": s1_ortholog_status["specific"]
        },
        f"{species_2_name}_ortholog_status": {
            "shared": s2_ortholog_status["shared"],
            "specific": s2_ortholog_status["specific"]
        }
    }

    
    #Step 3: Classify shared OCRs as promoters/enhancers
    logger.info("Step 3: Classifying shared OCRs")
    
    # Use species 1's shared regions and classify them using species 2 TSS
    shared_results = classifyOcrPromotersEnhancers(
        ocr_bed_path=s1_ortholog_status["shared"],
        tss_bed_path=config.get("species_2_tss_file_cleaned"),
        output_prefix=str(output_dir / "shared"),
        promoter_distance=promoter_distance,
        save_nearest=False
    )

    #Organizing results for logging
    shared_classified_results = {
        "shared_across_species": {
            "promoters": shared_results["promoters"],
            "enhancers": shared_results["enhancers"]
        }
    }


    #Step 4: Classify species-specific OCRs as promoters/enhancers
    logger.info("Step 4: Classifying species-specific OCRs")
    
    s1_specific_results = classifyOcrPromotersEnhancers(
        ocr_bed_path=s1_ortholog_status["specific"],
        tss_bed_path=config.get("species_1_tss_file_cleaned"),
        output_prefix=str(output_dir / f"{species_1_name}_specific"),
        promoter_distance=promoter_distance,
        save_nearest=False
    )
    
    s2_specific_results = classifyOcrPromotersEnhancers(
        ocr_bed_path=s2_ortholog_status["specific"],
        tss_bed_path=config.get("species_2_tss_file_cleaned"),
        output_prefix=str(output_dir / f"{species_2_name}_specific"),
        promoter_distance=promoter_distance,
        save_nearest=False
    )

    #Organizing results for logging
    specific_classified_results = {
        f"{species_1_name}_specific": {
            "promoters": s1_specific_results["promoters"],
            "enhancers": s1_specific_results["enhancers"]
        },
        f"{species_2_name}_specific": {
            "promoters": s2_specific_results["promoters"],
            "enhancers": s2_specific_results["enhancers"]
        }
    }
    

    #Print summary
    logger.info("Classification Complete")
    logger.info("Classification Results:")
    for species, files in classification_results.items():
        logger.info(f"  {species}:")
        for file_type, path in files.items():
            logger.info(f"    {file_type}: {path}")
    
    logger.info("Ortholog Status Classification:")
    for mapping, files in ortholog_results.items():
        logger.info(f"  {mapping}:")
        for file_type, path in files.items():
            logger.info(f"    {file_type}: {path}")
    
    logger.info("Shared Elements (Classified):")
    for category, files in shared_classified_results.items():
        logger.info(f"  {category}:")
        for file_type, path in files.items():
            logger.info(f"    {file_type}: {path}")
    
    logger.info("Species-Specific Elements (Classified):")
    for category, files in specific_classified_results.items():
        logger.info(f"  {category}:")
        for file_type, path in files.items():
            logger.info(f"    {file_type}: {path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Classify OCRs into promoters/enhancers and compare across species",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python3 classification.py --config config.processed.yaml
  python3 classification.py --config /path/to/config.processed.yaml 
  
Outputs: 10 BED files
    - {species}_all_promoters.bed
    - {species}_all_enhancers.bed
    - {species}_specific_promoters.bed
    - {species}_specific_enhancers.bed
    - shared_promoters.bed
    - shared_enhancers.bed
        """,
    )
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Path to YAML configuration file",
    )
    parser.add_argument(
        "--log-level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO)",
        )

    args = parser.parse_args()

    #Set log level
    logger.setLevel(getattr(logging, args.log_level))
    
    try:
        run_classification(args.config)
        logger.info("Pipeline finished successfully")
    except Exception as e:
        logger.error(f"Pipeline failed: {e}", exc_info=True)
        exit(1)
