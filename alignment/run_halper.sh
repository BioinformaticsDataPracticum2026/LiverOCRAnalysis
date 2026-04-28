#!/bin/bash
#SBATCH -J halper_map
#SBATCH -p RM-shared
#SBATCH -N 1
#SBATCH -n 2
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=2000M
#SBATCH -o logs/halper_%j.out
#SBATCH -e logs/halper_%j.err
#SBATCH --mail-type=END,FAIL

# usage:
# bash alignment/run_halper.sh human_peaks mouse_peaks hal_file output_dir

HUMAN_PEAKS=$1
MOUSE_PEAKS=$2
HAL_FILE=$3
OUTPUT_DIR=$4

# check inputs
if [ -z "$HUMAN_PEAKS" ] || [ -z "$MOUSE_PEAKS" ] || [ -z "$HAL_FILE" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: bash alignment/run_halper.sh <human_peaks> <mouse_peaks> <hal_file> <output_dir>"
    exit 1
fi

if [ ! -f "$HUMAN_PEAKS" ]; then
    echo "ERROR: human peaks file not found: $HUMAN_PEAKS"
    exit 1
fi

if [ ! -f "$MOUSE_PEAKS" ]; then
    echo "ERROR: mouse peaks file not found: $MOUSE_PEAKS"
    exit 1
fi

if [ ! -f "$HAL_FILE" ]; then
    echo "ERROR: HAL file not found: $HAL_FILE"
    exit 1
fi

# load environment
module load anaconda3
source activate hal

# create output folders
mkdir -p logs
mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/human_to_mouse"
mkdir -p "$OUTPUT_DIR/mouse_to_human"

# helper function for gz or plain bed files
prepare_bed () {
    INPUT_FILE=$1
    OUTPUT_FILE=$2

    if [[ "$INPUT_FILE" == *.gz ]]; then
        gunzip -c "$INPUT_FILE" > "$OUTPUT_FILE"
    else
        cp "$INPUT_FILE" "$OUTPUT_FILE"
    fi
}

# prepare peak files
echo "Preparing human peaks..."
prepare_bed "$HUMAN_PEAKS" "$OUTPUT_DIR/human_liver.narrowPeak"

echo "Preparing mouse peaks..."
prepare_bed "$MOUSE_PEAKS" "$OUTPUT_DIR/mouse_liver.narrowPeak"

echo "Starting HALPER mapping..."

# run human to mouse
bash alignment/halper.sh \
    -b "$OUTPUT_DIR/human_liver.narrowPeak" \
    -o "$OUTPUT_DIR/human_to_mouse" \
    -s "Human" \
    -t "Mouse" \
    -c "$HAL_FILE" &

# run mouse to human
bash alignment/halper.sh \
    -b "$OUTPUT_DIR/mouse_liver.narrowPeak" \
    -o "$OUTPUT_DIR/mouse_to_human" \
    -s "Mouse" \
    -t "Human" \
    -c "$HAL_FILE" &

wait

echo "HALPER mapping completed"