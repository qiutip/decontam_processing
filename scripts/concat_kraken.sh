#!/bin/bash

samples="$1"
file_output="$2"

awk '{ 
    file_path = $2;
    sample_id = $1;
    while ((getline line < file_path) > 0) {
        print sample_id, line;
    }
    close(file_path);
}' $samples > $file_output

