#!/bin/bash

path_out="$1"
kraken_dir="$2"

mkdir $path_out/output

python search_kraken.py -d $kraken_dir -o $path_out
kraken_samples="$path_out/output/kraken_samples.tsv"

echo "Concatenating all kraken.tsv (may take awhile)"
./concat_kraken.sh $kraken_samples $path_out

