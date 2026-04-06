from pybedtools import BedTool #Available in bioconda environment

#Human Liver:
#Load files
ocrHumanLiver = BedTool("ocrHumanLiver.bed")
tssHuman = BedTool("tssHuman.bed")

#Define ±2kb region in the promoter .bed file
promoterRegionsHuman = tssHuman.slop(genome="hg38.genome", b=2000)
promoterRegionsHuman.saveas("promoters2kbHuman.bed")

#Find promoters -> If overlap = Promoter
ocrPromotersHuman = ocrHumanLiver.intersect(promoterRegionsHuman, u=True)
ocrPromotersHuman.saveas("promotersHumanLiver.bed")

#Enhancers -> If no overlap = Enhancer
ocrEnhancersHuman = ocrHumanLiver.intersect(promoterRegionsHuman, v=True)
ocrEnhancersHuman.saveas("enhancersHumanLiver.bed")

#Find closest gene for each OCR
ocrWithGenesHuman = ocrHumanLiver.closest(tssHuman, d=True)
ocrWithGenesHuman.saveas("ocrNearestTSSHuman.bed")



#Mouse Liver:
#Load files
ocrMouseLiver = BedTool("ocrMouseLiver.bed")
tssMouse = BedTool("tssMouse.bed")

#Define ±2kb region in the promoter .bed file
promoterRegionsMouse = tssMouse.slop(genome="mm10.genome", b=2000)
promoterRegionsMouse.saveas("promoters2kbMouse.bed")

#Find promoters -> If overlap = Promoter
ocrPromotersMouse = ocrMouseLiver.intersect(promoterRegionsMouse, u=True)
ocrPromotersMouse.saveas("promotersMouseLiver.bed")

#Enhancers -> If no overlap = Enhancer
ocrEnhancersMouse = ocrMouseLiver.intersect(promoterRegionsMouse, v=True)
ocrEnhancersMouse.saveas("enhancersMouseLiver.bed")

#Find closest gene for each OCR
ocrWithGenesMouse = ocrMouseLiver.closest(tssMouse, d=True)
ocrWithGenesMouse.saveas("ocrNearestTSSMouse.bed")



#Probably better to create a wrapper function:
"""
classifyOcrPromotersEnhancers() --> Classify OCRs as promoters or enhancers and assign nearest TSS.
Paramaeters:
- ocrBedPath (str): Path to the OCR BED file (peaks).
- tssBedPath (str): Path to the TSS BED file (1bp per TSS).
- genomeFile (str): Genome file (chromosome sizes) for pybedtools slop.
- outputPrefix (str): Prefix for output files. Outputs will be saved as:
    - {outputPrefix}_promoters.bed
    - {outputPrefix}_enhancers.bed
    - {outputPrefix}_nearestTSS.bed
 """
def classifyOcrPromotersEnhancers(
    ocrBedPath: str,
    tssBedPath: str,
    genomeFile: str,
    outputPrefix: str
):
    
    # Load BED files
    ocr = BedTool(ocrBedPath)
    tss = BedTool(tssBedPath)
    
    # Expand TSS to ±2kb promoter regions
    promoter_regions = tss.slop(genome=genomeFile, b=2000)
    promoter_regions.saveas(f"{outputPrefix}_promoterRegions2kb.bed")
    
    # Classify promoters (overlap)
    ocrPromoters = ocr.intersect(promoter_regions, u=True)
    ocrPromoters.saveas(f"{outputPrefix}_promoters.bed")
    
    # Classify enhancers (no overlap)
    ocrEnhancers = ocr.intersect(promoter_regions, v=True)
    ocrEnhancers.saveas(f"{outputPrefix}_enhancers.bed")
    
    # Assign nearest TSS
    ocrWithGenes = ocr.closest(tss, d=True)
    ocrWithGenes.saveas(f"{outputPrefix}_nearestTSS.bed")
    
    print(f"Finished processing {ocrBedPath}.")
    print(f"Promoters: {outputPrefix}_promoters.bed")
    print(f"Enhancers: {outputPrefix}_enhancers.bed")
    print(f"Nearest TSS: {outputPrefix}_nearestTSS.bed\n")
    
    return ocrPromoters, ocrEnhancers, ocrWithGenes

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Classify OCRs into promoters/enhancers")
    parser.add_argument("--config", required=True)

    args = parser.parse_args()

    run_classification(Path(args.config))
