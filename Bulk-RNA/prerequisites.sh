#!/bin/bash

base_dir="command_inputs"

# Function to check the success of a command
check_success() {
  if [ $? -ne 0 ]; then
    echo "Command '$1' failed."
    exit 1
  fi
}

# Create directories
mkdir -p "$base_dir"/{Aligner-index,RNA-GTF,adapters,ACTB,tidy}
check_success "mkdir"

# Change to the base directory
cd "$base_dir" || exit 1

# Check if the RNA GTF file exists
if [ -f "./RNA-GTF/hs1.ncbiRefSeq.gtf" ]; then
  echo "RNA GTF file exists; no downloading will be done."
else
  # Download and extract the RNA GTF file
  wget -c "https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/genes/hs1.ncbiRefSeq.gtf.gz" -O ./RNA-GTF/hs1.ncbiRefSeq.gtf.gz
  check_success "wget"
  
  gunzip ./RNA-GTF/hs1.ncbiRefSeq.gtf.gz
  check_success "gunzip"
fi

# Check if the ACTB BED file exists
if [ -f "./ACTB/ACTB.bed" ]; then
  echo "BED file exists; no downloading will be done."
else
  echo -e "chr7\t5527147\t5530601\tNM_001101.5\t0\t-\t5527747\t5529657\t0\t6\t744,182,439,240,129,78,\t0,856,1133,2013,2387,3376," > ./ACTB/ACTB.bed
fi

# Check if the HISAT2 index files exist
if [ -f "./Aligner-index/genome.1.ht2" ]; then
  echo "Index files exist; no downloading will be done."
else
  # Download and extract the HISAT2 index files for hg38
  wget -c "https://genome-idx.s3.amazonaws.com/hisat/hg38_genome.tar.gz" -O ./Aligner-index/hg38_genome.tar.gz
  check_success "wget"
  
  tar -xvzf ./Aligner-index/hg38_genome.tar.gz -C ./Aligner-index/
  check_success "tar"
fi

# Check if the tidy script exists
if [ -f "./tidy/stringtie_expression_matrix.pl" ]; then
  echo "Tidy script exists; no downloading will be done."
else
  wget -c "https://raw.githubusercontent.com/griffithlab/rnabio.org/master/assets/scripts/stringtie_expression_matrix.pl" -O ./tidy/stringtie_expression_matrix.pl
  check_success "wget"
  
  chmod +x ./tidy/stringtie_expression_matrix.pl
fi

# Check if adapter sequences exist
if [ -f "./adapters/TruSeq2-PE.fa" ]; then
  echo "Adapters exist; no downloading will be done."
else
  wget -c "https://raw.githubusercontent.com/timflutre/trimmomatic/master/adapters/TruSeq2-PE.fa" -O ./adapters/TruSeq2-PE.fa
  check_success "wget"
fi


