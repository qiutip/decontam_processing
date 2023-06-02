#!/bin/bash

tax_file="$1"
kraken_file="$2"
output_file="$3"

awk -v tax=$tax_file 'BEGIN {
    comm="cat "tax;
    while(comm|getline) {
        conditions[$1]=1;
    }
}{
    {
        if($4 in conditions) {
            interest=1;
        } else {
            interest=0;
        }
    } if (!interest==1) {
        print $1,"@"$3;}
}' $kraken_file > $output_file
