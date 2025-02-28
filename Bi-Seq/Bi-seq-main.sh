#!/bin/bash

base_dir="GSE1011"
srr_file="$HOME/Downloads/command-inputs/SRR_Acc_List.txt"
adapter="$HOME/Downloads/command-inputs/adapter/"
refrenence="$HOME/Downloads/command-inputs/reference-index/"
thread_numbers="4"


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
mkdir -p "$base_dir"/{fastq/{trimmed,untrimmed},fastqc/{trimmed,untrimmed},alignment/{raw,sorted},methylation-profile}
check_success "mkdir -p"

while IFS= read -r line; do
    # Store each line (SRA ID) into the SraAcc variable
    SraAcc="$line"

  for i in ./fastq/untrimmed/"${SraAcc}"_1.fastq.gz; do
      # Define output file names based on the input file
      base=$(basename "$i" _1.fastq.gz)
      paired1="./fastq/trimmed/paired/${base}_1_paired.fastq.gz"
      unpaired1="./fastq/trimmed/unpaired/${base}_1_unpaired.fastq.gz"
      paired2="./fastq/trimmed/paired/${base}_2_paired.fastq.gz"
      unpaired2="./fastq/trimmed/unpaired/${base}_2_unpaired.fastq.gz"

      # Define the second input file
      second_input="./fastq/untrimmed/${base}_2.fastq.gz"

      # Check if both input files exist
      if [[ ! -f "$second_input" ]]; then
          echo "Second input file $second_input does not exist. Skipping Trimmomatic for $SraAcc."
          continue  # Skip to the next iteration if the second input file does not exist
      fi

    # Check if output files already exist
      if [[ -f "$paired1" && -f "$unpaired1" && -f "$paired2" && -f "$unpaired2" ]]; then
        echo "Output files already exist for $SraAcc. Skipping Trimmomatic."
        continue  # Skip to the next iteration if output files exist
      fi
      
      # Run Trimmomatic if input and output files are valid
      trimmomatic PE "$i" "$second_input" \
        "$paired1" \
        "$unpaired1" \
        "$paired2" \
        "$unpaired2" \
        ILLUMINACLIP:"$adapter":2:30:10 \
        -threads "$thread_numbers"

      check_success "trimmomatic"
  done
   
    #fastqc
    for i in ./fastq/trimmed/paired/"${SraAcc}"_*.fastq.gz; do
        fastqc "$i" --threads "$thread_numbers" --outdir ./fastqc/trimmed/
        check_success "fastqc"
    done
    
    
   for i in ./fastq/trimmed/paired/*_1_paired.fastq.gz; do
    # Extract the base name without directory and suffix
      base=$(basename "$i" _1_paired.fastq.gz)

      # Define expected Bowtie2 output file name
      output="./alignment/raw/${base}.bam"

      # Check if the SAM file already exists
      if [[ -f "$sam_output" ]]; then
        echo "alignment output already exists for $base. Skipping alignment."
        continue  # Skip to the next iteration if the SAM file exists
      fi
      
      echo "performing alignment for $base"
      # Run Bowtie2 if SAM file does not exist
      bismark \
        "$reference" \
        -1 "$i" \
        -2 "./fastq/trimmed/paired/${base}_2_paired.fastq.gz" \
        -o "$output" \
        --parallel "$thread_numbers"
      echo "finished alignment for $base \n" 
      check_success "bismark"
  done
  
  for i in ./alignment/raw/*.bam; do
    # Extract the base name without directory and suffix
     base=$(basename "$i" .bam)

     mkdir "./methylation-profile/$base"

     bismark_methylation_extractor \
     --cytosine_report \
     $i \
     --genome_folder "$refrenence"
     -o "./methylation-profile/$base/"
     --parallel "$thread_number"
      
     check_success "bowtie2"
  done
done < "$srr_file" | tee "$HOME/Downloads/$base_dir/README.md" 2>&1

