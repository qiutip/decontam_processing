# decontam_processing
filters out likely contanimated sequences using decontam and taxonomy data (e.g. Kraken2)

Step 1.
use scripts/process_kraken_total.sh path/to/directory path/to/directory/kraken_samples.txt

Step 2.
Use scripts/process.sh path/to/directory path/to/taxaonomy/merged/bracken_species.tsv path/to/merged/bracken_metadata.tsv path/to/kraken_total.tsv path/to/samples_fastq.tsv
