#!/bin/bash

# Array of inputs
inputs=("/athena/masonlab/scratch/users/jaq4005/decontam_processing/test/samples_fastq.txt" "/athena/masonlab/scratch/users/jaq4005/stuckonu/notebook/decontam_data/test_concat/output/read_id.filter.txt" "/athena/masonlab/scratch/users/jaq4005/decontam_processing/test/parallel")

#path_samples=$1
#path_read_id=$2
#path_output=$3

# Export the function to make it accessible by parallel
#export -f extract_fastq_kraken.sh

# Run the script in parallel using GNU Parallel
#./extract_fastq_kraken.sh $path_samnples $path_read_id $path_output | parallel
parallel ./extract_fastq_kraken.sh ::: "${inputs[0]}" ::: "${inputs[1]}" ::: "${inputs[2]}" 
#parallel -j 4 ./extract_fastq_kraken.sh {} {} {} ::: $path_samples ::: $path_read_id ::: $path_output
#parallel -j 4 ./extract_fastq_kraken.sh "${inputs[0]}" "${inputs[1]}" "${inputs[2]}"
