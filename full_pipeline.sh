#!/bin/bash
#SBATCH -J liver_ocr_pipeline
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH -o logs/pipeline_%j.out
#SBATCH -e logs/pipeline_%j.err
#SBATCH --mail-type=END,FAIL

# create logs folder
mkdir -p logs

# load HOMER for motif step
module load homer

# -------- alignment --------
python main.py alignment \
  --human-peaks data/atac_seq_peaks/humanIDRConservedPeaks.gz \
  --mouse-peaks data/atac_seq_peaks/mouseIDRConservedPeaks.gz \
  --hal-file /path/to/10plusway-master.hal \
  --outdir results/alignment_results \
  --local

# -------- classification --------
python main.py classification \
  --config classification/sample_config.yaml

# -------- motif --------
python main.py motif \
  --genome hg38 \
  --bed-dir results/classification_results \
  --beds human_all_promoters.bed human_all_enhancers.bed shared_promoters.bed shared_enhancers.bed human_specific_promoters.bed human_specific_enhancers.bed

# -------- enrichment --------
python -m enrichment_analysis.run_great