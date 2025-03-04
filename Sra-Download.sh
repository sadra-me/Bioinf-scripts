#!/bin/bash

# Define variables
base_dir="GSE1010"
srr_file="$HOME/Downloads/command_inputs/SRR_Acc_List.txt"  # Updated path
thread_numbers="4"

# Function to check if a command was successful
check_success() {
    if [ $? -eq 0 ]; then
        echo "Command '$1' succeeded."
    else
        echo "Command '$1' failed."
        exit 1  # Exit the script if a critical command fails
    fi
}

# Create the base directory and subdirectories
mkdir -p "$base_dir"/{fastq/{trimmed,untrimmed},fastqc/{trimmed/{paired,unpaired},untrimmed}}
check_success "mkdir -p"

# Check if the base directory was created successfully
if [ ! -d "$base_dir" ]; then
    echo "Error: Base directory '$base_dir' was not created."
    exit 1
fi

# Change the current directory to the base directory once
cd "$base_dir" || { echo "Could not change directory to '$base_dir'"; exit 1; }

# Read the file line by line
while IFS= read -r line; do
    # Store each line (SRA ID) into the SraAcc variable
    SraAcc="$line"

    # Check if output files already exist (either .fastq or .fastq.gz)
    if [[ -f "./fastq/untrimmed/${SraAcc}_1.fastq" && -f "./fastq/untrimmed/${SraAcc}_2.fastq" ]] || \
       [[ -f "./fastq/untrimmed/${SraAcc}_1.fastq.gz" && -f "./fastq/untrimmed/${SraAcc}_2.fastq.gz" ]]; then
        echo "Output files for '$SraAcc' already exist. Skipping prefetch and fastq-dump."
        continue
    fi

    # Use prefetch command with the SraAcc variable
    prefetch "$SraAcc" --output-directory ./fastq/untrimmed/
    check_success "prefetch"
    
    vdb-validate "./fastq/untrimmed/$SraAcc/$SraAcc.sra" | tee "./fastq/validation.log"
    check_success "vdb-validate"

    # Perform fastq-dump with --split-files
    fasterq-dump -e "$thread_numbers" --split-files "$SraAcc" -O ./fastq/untrimmed/
    check_success "fasterq-dump"

    # Perform gzip on the resulting fastq files using pigz, only if they are not already gzipped
    if [[ ! -f "./fastq/untrimmed/${SraAcc}_1.fastq.gz" || ! -f "./fastq/untrimmed/${SraAcc}_2.fastq.gz" ]]; then
        pigz -p "$thread_numbers" ./fastq/untrimmed/"${SraAcc}"_*.fastq
        check_success "pigz"
    fi

    # Run FastQC on the gzipped fastq files
    for i in ./fastq/untrimmed/"${SraAcc}"_1.fastq.gz; do
        fastqc "$i" --threads "$thread_numbers" --outdir ./fastqc/untrimmed/
        fastqc "${i/_1.fastq/_2.fastq}" --threads "$thread_numbers" --outdir ./fastqc/untrimmed/
        check_success "fastqc"
    done

done < "$srr_file" | tee "$HOME/Downloads/$base_dir/fastq/README.md" 2>&1

