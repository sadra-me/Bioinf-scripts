#!/bin/bash

base_dir="GSE1010"
srr_file="$HOME/Downloads/command_inputs/SRR_Acc_List.txt"
adapter="$HOME/Downloads/command_inputs/adapters/TruSeq2-PE.fa"
reference="$HOME/Downloads/command_inputs/Aligner-index/hg38/genome"
GeneBodyCov_ACTB="$HOME/Downloads/command_inputs/ACTB/ACTB.bed"
RNA_GTF="$HOME/Downloads/command_inputs/RNA-GTF/hs1.ncbiRefSeq.gtf"
Tidy_loc="$HOME/Downloads/command_inputs/tidy/stringtie_expression_matrix.pl"
thread_numbers="4"

# Check if base directory exists
if [ ! -d "$base_dir" ]; then
  echo "Error: Base directory '$base_dir' not found."
  exit 1
fi

# Check if SRR accession file exists and is not empty
if [[ ! -s "$srr_file" ]]; then
    echo "Error: SRR accession file '$srr_file' is missing or empty."
    exit 1
fi

cd "$base_dir"

# Function to check command success and handle errors
check_success() {
  if [ $? -ne 0 ]; then
    echo "Command '$1' failed."
    exit 1
  fi
}

# Create necessary directories
mkdir -p alignment/log QC-results/{untrimmed,RSeQC} expression/StringTie/
check_success "mkdir"

# Process each SRR accession in the list
while IFS= read -r line; do
    SraAcc="$line"


    # Hisat2 Alignment Loop (Alignment)
    for i in ./fastq/untrimmed/*_1.fastq.gz; do
        base=$(basename "$i" _1.fastq.gz)
        output="./alignment/${base}.sam"
	read_1="./fastq/untrimmed/${base}_1.fastq.gz"
	read_2=".//fastq/untrimmed/${base}_2.fastq.gz"
        if [[ -f "$output" || -f "./alignment/${base}-sorted.bam" ]]; then
            echo "Alignment output already exists for $base. Skipping alignment."
            continue
        fi

        hisat2 \
            -x "$reference" \
            -1 "$read_1" \
            -2 "$read_2" \
            -S "$output" \
            --summary-file "./alignment/log/${base}_summary.txt" \
            -p "$thread_numbers"
        check_success "hisat2"
    done

    # SAMtools Sorting Loop (Sorting BAM Files)
    for i in ./alignment/*.sam; do
        base=$(basename "$i" .sam)
        output="./alignment/${base}-sorted.bam"

        if [[ -f "$output" ]]; then
            echo "${base}-sorted.bam already exists. Skipping sorting."
            rm ./alignment/"${base}.sam"
            continue
        fi

        samtools sort -@ "$thread_numbers" "$i" -o "$output"
        check_success "samtools sort"
    done

    # StringTie Expression Analysis Loop (Quantification)
    for i in ./alignment/*-sorted.bam; do
        base=$(basename "$i" -sorted.bam)
        
        mkdir -p ./expression/StringTie/"$base"
        out_dir="./expression/StringTie/$base/"
        
        stringtie \
            $i \
            --rf \
            -p "$thread_numbers" \
            -G $RNA_GTF \
            -e -B \
            -o "${out_dir}/transcripts.gtf" \
            -A "${out_dir}/gene_abundance.tsv"
        
        check_success "stringtie"
    done

done < "$srr_file"

# Generate Expression Matrices with Tidy Script (Post-Processing)
directories=$(find ./expression/StringTie/ -mindepth 1 -maxdepth 1 -type d | tr '\n' ',' | sed 's/,$//')

"$Tidy_loc" \
    --expression_metric=TPM \
    --result_dirs="$directories" \
    --transcript_matrix_file=./expression/transcript_tpm_all_samples.tsv \
    --gene_matrix_file=./expression/gene_tpm_all_samples.tsv

check_success "TPM matrix generation"

"$Tidy_loc" \
    --expression_metric=FPKM \
    --result_dirs="$directories" \
    --transcript_matrix_file=./expression/transcript_fpkm_all_samples.tsv \
    --gene_matrix_file=./expression/gene_fpkm_all_samples.tsv

check_success "FPKM matrix generation"

"$Tidy_loc" \
    --expression_metric=Coverage \
    --result_dirs="$directories" \
    --transcript_matrix_file=./expression/transcript_coverage_all_samples.tsv \
    --gene_matrix_file=./expression/gene_coverage_all_samples.tsv

check_success "Coverage matrix generation"

