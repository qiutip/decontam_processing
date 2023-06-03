#!/bin/bash

path_output=$1
path_bracken_combined=$2
path_bracken_meta=$3
path_kraken_total=$4
path_samples_fastq=$5

mkdir $path_output/output
mkdir $path_output/output/decontam

echo "Creating Tax_ids to filter"
Rscript ./decontam_processing.R $path_output $path_bracken_combined $path_bracken_meta

echo "Extracting Read IDs to filter"
./extract_reads_kraken.sh $path_output"/output/tax_list.txt" $path_kraken_total $path_output"/output/read_id.filter.txt" 

echo "Removing likely contanimated sequences from FastQs"
./extract_fastq_kraken.sh $path_samples_fastq $path_output"/output/read_id.filter.txt" $path_output
