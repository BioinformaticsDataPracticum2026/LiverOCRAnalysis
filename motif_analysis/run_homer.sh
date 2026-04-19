#!/bin/bash
#SBATCH -J enrichment_all
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH -n 8  
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=smakkar@andrew.cmu.edu

export PATH=$PATH:/jet/home/smakkar/enrichment/homer/bin/

#directories
INPUT_DIR="/jet/home/smakkar/enrichment/raw_results"
BASE_OUT="/jet/home/smakkar/enrichment/results"

mkdir -p $BASE_OUT/human_promoters
mkdir -p $BASE_OUT/human_enhancers

echo "Starting all HOMER analyses concurrently"

# Human Promoters (hg38)
findMotifsGenome.pl $INPUT_DIR/human_promoters.bed hg38 $BASE_OUT/human_promoters -size 200 -mask -p 2 &

# Human Enhancers (hg38)
findMotifsGenome.pl $INPUT_DIR/human_enhancers.bed hg38 $BASE_OUT/human_enhancers -size 200 -mask -p 2 &


wait

echo "All motif analyses complete!"