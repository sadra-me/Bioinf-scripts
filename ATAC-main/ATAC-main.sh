#!/bin/bash

# Set -e: Exit immediately if a command exits with a non-zero status
set -e

# Define variables
base_dir="GSE1010"
srr_file="$HOME/Downloads/command-inputs/SRR_Acc_List.txt"
adapter="$HOME/Downloads/command-inputs/adapter/NexteraPE-PE.fa"
reference="$HOME/Downloads/command-inputs/reference-index/GRCh38_noalt_as/GRCh38_noalt_as"
blacklist="$HOME/Downloads/command-inputs/blacklist/ENCFF356LFX.bed"
chrom_size="$HOME/Downloads/command-inputs/chrom-size/hg38.chrom.sizes"
thread_numbers="4"

# Check if the base directory exists
if [ ! -d "$base_dir" ]; then
  echo "Error: Base directory '$base_dir' not found."
  exit 1
fi

# Change to the base directory
cd "$base_dir"

# Define check_success function
check_success() {
  if [ $? -ne 0 ]; then
    echo "Command '$1' failed."
    exit 1
  fi
}

# Ensure output directories exist
mkdir -p ./peaks/blacklisted
mkdir -p ./bigwig

# BedClip and BigWig Conversion
for i in ./peaks/blacklisted/*-blacklist-Removed.bdg; do
  # Extract the base name
  base=$(basename "$i" "-blacklist-Removed.bdg")

  # Define the output file names
  clipped_output="./peaks/blacklisted/${base}-blacklist-Removed-clipped.bdg"
  bigwig_output="./bigwig/${base}.bw"

  # Check if both the clipped file and the BigWig file already exist
  if [[ -f "$clipped_output" && -f "$bigwig_output" ]]; then
    echo "Both clipped and BigWig files exist for ${base}. Skipping."
    continue  # Skip to the next iteration
  fi

  # Perform BedClip if the clipped file doesn't exist
  if [[ ! -f "$clipped_output" ]]; then
    echo "Performing bedClip for $i"
    bedClip "$i" "$chrom_size" "$clipped_output"
    check_success "bedClip"
  fi

  # Convert to BigWig if the BigWig file doesn't exist
  if [[ ! -f "$bigwig_output" ]]; then
    echo "Performing BigWig conversion for $i"

    # Define temporary output path
    temp_bdg="./bigwig/tmp.bdg"

    # Sort the bedGraph file and save to temporary file
    sort -k1,1 -k2,2n "$clipped_output" > "$temp_bdg"

    # Convert bedGraph to BigWig
    bedGraphToBigWig "$temp_bdg" "$chrom_size" "$bigwig_output"

    # Remove the temporary file
    rm "$temp_bdg"

    check_success "bedGraphToBigWig"
  fi
done

