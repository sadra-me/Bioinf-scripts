#!/bin/bash

base_dir="GSE1011"
srr_file="$HOME/Downloads/command-inputs/SRR_Acc_List.txt"
adapter="$HOME/Downloads/command-inputs/adapter/"
reference="$HOME/Downloads/command-inputs/reference-index/"
GeneBodyCov_ACTB="$HOME/Downloads/command-inputs/genebody/ACTB.bed"
RNA_GTF="$HOME/Downloads/command-inputs/gtf/reference.gtf"
Tidy_loc="$HOME/Downloads/command-inputs/scripts/stringtie_expression_matrix.pl"
thread_numbers="4"

if [ ! -d "$base_dir" ]; then
  echo "Error: Base directory '$base_dir' not found."
  exit 1
fi

cd "$base_dir"

check_success() {
  if [ $? -ne 0 ]; then
    echo "Command '$1' failed."
    exit 1
  fi
}

mkdir -p fastq/{trimmed,untrimmed} fastqc/{trimmed,untrimmed} alignment/{raw,sorted} \
methylation-profile QC-results/{trimmed,untrimmed,RSeQC} expression/StringTie/
check_success "mkdir"

while IFS= read -r line; do
    SraAcc="$line"

    # Trimmomatic Loop
    for i in ./fastq/untrimmed/"${SraAcc}"_1.fastq.gz; do
        base=$(basename "$i" _1.fastq.gz)
        paired1="./fastq/trimmed/${base}_1_paired.fastq.gz"
        unpaired1="./fastq/trimmed/${base}_1_unpaired.fastq.gz"
        paired2="./fastq/trimmed/${base}_2_paired.fastq.gz"
        unpaired2="./fastq/trimmed/${base}_2_unpaired.fastq.gz"
        second_input="./fastq/untrimmed/${base}_2.fastq.gz"

        if [[ ! -f "$i" || ! -f "$second_input" ]]; then
            echo "Input files $i or $second_input do not exist. Skipping Trimmomatic for $SraAcc."
            continue
        fi

        if [[ -f "$paired1" && -f "$unpaired1" && -f "$paired2" && -f "$unpaired2" ]]; then
            echo "Output files already exist for $SraAcc. Skipping Trimmomatic."
            continue
        fi

        trimmomatic PE "$i" "$second_input" \
            "$paired1" "$unpaired1" \
            "$paired2" "$unpaired2" \
            ILLUMINACLIP:"$adapter":2:30:10 \
            -threads "$thread_numbers"
        check_success "trimmomatic"
    done

    # FastQC and MultiQC Loop
    for i in ./fastq/trimmed/"${SraAcc}"_*.fastq.gz; do
        fastqc "$i" --threads "$thread_numbers" --outdir ./fastqc/trimmed/
        check_success "fastqc"
    done

    multiqc ./fastqc/trimmed/ -o ./QC-results/trimmed/
    check_success "multiqc"

    # Hisat2 Alignment Loop
    for i in ./fastq/trimmed/*_1_paired.fastq.gz; do
        base=$(basename "$i" _1_paired.fastq.gz)
        output="./alignment/raw/${base}.sam"

        if [[ -f "$output" ]]; then
            echo "Alignment output already exists for $base. Skipping alignment."
            continue
        fi

        hisat2 \
            -x "$reference" \
            -1 "$i" \
            -2 "./fastq/trimmed/${base}_2_paired.fastq.gz" \
            -S "$output" \
            --summary-file "./alignment/log/${base}_summary.txt" \
            -p "$thread_numbers"
        check_success "hisat2"
    done

    # SAMtools Sorting Loop
    for i in ./alignment/raw/*.sam; do
        base=$(basename "$i" .sam)
        output="./alignment/sorted/${base}-sorted.bam"

        if [[ -f "$output" ]]; then
            echo "${base}-sorted.bam already exists. Skipping sorting."
            continue
        fi

        samtools sort -@ "$thread_numbers" "$i" -o "$output"
        check_success "samtools sort"
    done

    # StringTie Expression Analysis Loop
    for i in ./alignment/sorted/*-sorted.bam; do
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

# Generate Expression Matrices with Tidy Script
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

echo "Pipeline completed successfully!"

