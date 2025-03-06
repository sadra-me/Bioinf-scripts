#!/bin/bash
set -euo pipefail

# Define variables
base_dir="GSE1010"
srr_file="$HOME/Downloads/command_inputs/SRR_Acc_List.txt"
thread_numbers="4"

# Function to check if a command was successful
check_success() {
    if [ $? -eq 0 ]; then
        echo "Command '$1' succeeded."
    else
        echo "Command '$1' failed."
        exit 1
    fi
}

# Create directories
mkdir -p "$base_dir"/{fastq/{trimmed,untrimmed},fastqc/{trimmed/{paired,unpaired},untrimmed}}


cd "$base_dir"


while IFS= read -r SraAcc; do
    # Download SRA file with prefetch
    prefetch "$SraAcc" --output-directory ./fastq/untrimmed/
    check_success "prefetch"
done < "$srr_file"


while IFS= read -r SraAcc; do
    # Check if output FASTQ files already exist
    if [[ -f "./fastq/untrimmed/${SraAcc}_1.fastq.gz" && -f "./fastq/untrimmed/${SraAcc}_2.fastq.gz" ]]; then
        echo "Output files for '$SraAcc' already exist. Skipping processing."
        continue
    fi

    # Validate downloaded data
    vdb-validate "./fastq/untrimmed/$SraAcc/$SraAcc.sra" | tee -a "./fastq/validation.log"
    check_success "vdb-validate"

    # Convert to FASTQ format with fasterq-dump
    fasterq-dump -e "$thread_numbers" --split-files "$SraAcc" -O ./fastq/untrimmed/
    check_success "fasterq-dump"

    # Compress FASTQ files using pigz
    pigz -p "$thread_numbers" ./fastq/untrimmed/"${SraAcc}"_*.fastq
    check_success "pigz"

    # Run FastQC on compressed FASTQ files
    fastqc ./fastq/untrimmed/"${SraAcc}"_1.fastq.gz --threads "$thread_numbers" --outdir ./fastqc/untrimmed/
    fastqc ./fastq/untrimmed/"${SraAcc}"_2.fastq.gz --threads "$thread_numbers" --outdir ./fastqc/untrimmed/
    check_success "fastqc"
done < "$srr_file"

echo "Processing phase completed."

