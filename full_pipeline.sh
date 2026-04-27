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

# create logs folder if it does not exist
mkdir -p logs

# run alignment (uses SLURM internally)
python main.py alignment

# run classification (includes preprocessing)
python main.py classification \
  --config classification/sample_config.yaml

# run motif analysis
python main.py motif \
  --genome hg38

# run enrichment analysis
python main.py enrichment