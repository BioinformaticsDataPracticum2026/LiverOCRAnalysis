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
python3 classification.py --config config.yaml