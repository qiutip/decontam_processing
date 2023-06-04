# decontam_processing
filters out likely contanimated sequences using decontam and taxonomy data (e.g. Kraken2)

Step 1.
use scripts/process_kraken_total.sh path/to/directory path/to/directory/kraken_samples.txt

Step 2.
Use scripts/process.sh path/to/directory path/to/taxaonomy/merged/bracken_species.tsv path/to/merged/bracken_metadata.tsv path/to/kraken_total.tsv path/to/samples_fastq.tsv


#Executing
## 1. process_kraken_total.sh
        Concatenates all of the kraken.tsv files with associated read_ids
                Creates output/kraken_samples.txt
                Creates output/kraken_total.txt
## 2. process.sh
        Identifies likely contanimated taxa
                Creates output/tax_list.txt
        Identifies read_ids associasted with cotanimated taxa
                Creates output/read_id.filter.txt
        Filters out FASTQ files
                Creates output/decontam/{sample_name}.R1.fastq.gz
                Creates output/decontam/{sample_name}.R2.fastq.gz

### Example
./process_kraken_total.sh $test_dir $kraken_dir

$kraken_dir=/athena/masonlab/scratch/users/jaq4005/stuckonu/bioforward/short-read-taxonomy/controls/short-read-taxonomy/2_kraken2/raw_kraken: Is a directory
$test_dir=/athena/masonlab/scratch/users/jaq4005/stuckonu/notebook/decontam_data/test_concat: Is a directory

### Example
./process.sh $test_dir $bracken_species_file $bracken_meta $kraken_total $samples_fastq 

$bracken_species_file=/athena/masonlab/scratch/users/jaq4005/stuckonu/bioforward/short-read-taxonomy/controls/short-read-taxonomy/2_kraken2/merged/bracken_species.tsv
$bracken_meta=/athena/masonlab/scratch/users/jaq4005/stuckonu/notebook/decontam_data/output/df_bs.meta.tsv
$kraken_total=/athena/masonlab/scratch/users/jaq4005/stuckonu/notebook/decontam_data/test_concat/output/kraken_total.tsv
$samples_fastq=/athena/masonlab/scratch/users/jaq4005/stuckonu/notebook/decontam_data/test_concat/samples.low.txt
