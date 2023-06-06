# decontam_processing

## Step 1.
Make a meta_table of the samples dictated as "Control" or "True" with column_name = "Sample_or_Control"
    Look at input for examples.

## Step 2.
1. use Rscript scripts/decontam_process.v3.R $path_output $path_otu $path_meta $dtype $threshold_interval $workers
    - dtype = "kraken" or "xtree"
        - Will implement metaplhan
    - threshold_interval = interval1,interval2
        - input as a tuple, seperator = ","
        - interval1 and interval2 are starting and ending points of seq(interval1,interval2, by = 5)

## Outputs
| Files/Figures |
| --- |
| generates tsv file with new relative abundances based on filtering out specific taxa across different taxaonomic ranks |
| generates tsv file with taxa filtered out at specific thresholds |
| generates tsv file of threshold vs contaminated taxa |
| produces a figure of threshold vs contaminated taxa |


## EXAMPLE
    Rscript scripts/decontam_process.v3.R path/to/output_directory path/to/input/xtree/GTDB_.1_.05_c_metagenomics_ra.tsv path/to/input/xtree/df.meta.tsv xtree .1,.5 5
    Rscript scripts/decontam_processing.v3.R /athena/masonlab/scratch/users/jaq4005/decontam_processing input/xtree/GTDB_.1_.05_c_metagenomics_ra.tsv input/xtree/df.meta.tsv xtree .1,.5 5

