#!/bin/bash

##create for loop that does both r1 and r2 based on samples.csv file

echo "Removing Contanimated Sequences"
file_name="$1"
read_id_file="$2"
path_out="$3"

sample_list=$(awk '{ print $1 }' $file_name)
read1_list=$(awk '{ print $2 }' $file_name)
read2_list=$(awk '{ print $3 }' $file_name)

samples=()
for element in $sample_list; do
    samples+=("$element")
done

read1=()
for element in $read1_list; do
    read1+=("$element")
done

read2=()
for element in $read2_list; do
    read2+=("$element")
done

for index in "${!samples[@]}"; do
    
    sample_val="${samples[$index]}" 
    read1_val="${read1[$index]}"
    read2_val="${read2[$index]}"
    dir_1="$path_out/output/decontam/"$sample_val"_1.fastq.gz"
    dir_2="$path_out/output/decontam/"$sample_val"_2.fastq.gz"
    
    echo "Filtering "$sample_val" Read 1"
    zcat $read1_val |  awk -v id=$read_id_file -v sample_id=$sample_val 'BEGIN {
        comm="cat "id;
        while(comm|getline) {
            if($1 == sample_id) {
                ids[$2]=1
                }
            }
        }{
            n++; 
                if (n % 4==1) {
                    if ($1 in ids) {
                        interest=1;
                    } 
                    else {
                        interest=0;}
                    }
                    if (interest==1) {
                        print;
                    }
        }' | gzip > "$dir_1"

    echo "Filtering "$sample_val" Read 2"    
    zcat $read2_val |  awk -v id=$read_id_file -v sample_id=$sample_val 'BEGIN {
        comm="cat "id;
        while(comm|getline) {
            if($1 == sample_id) {
                ids[$2]=1
                }
            }
        }{
            n++; 
                if (n % 4==1) {
                    if ($1 in ids) {
                        interest=1;
                    } 
                    else {
                        interest=0;}
                    }
                    if (interest==1) {
                        print;
                    }
        }' | gzip > "$dir_2" 
done
