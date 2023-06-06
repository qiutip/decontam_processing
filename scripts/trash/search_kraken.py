import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='Extract kraken.tsv files from CAMP module.')
parser.add_argument('-d', '--directory', help='Path to directory for samples before kraken.tsv files', required=True)
parser.add_argument('-o', '--output_dir', help='Path to out', required=True)
args = parser.parse_args()
directory = args.directory
output_path = args.output_dir

file_list = []
sample_list = []
for filename in os.listdir(directory):
    d = os.path.join(directory, filename)
    if os.path.isdir(d):
        for filename_2 in os.listdir(d):
            f = os.path.join(d, filename_2)
            if f.split("/")[-1] == "kraken.tsv":
                file_list.append(f)
                sample_list.append(f.split("/")[-2])

df_file_list = pd.DataFrame(zip(sample_list,file_list), columns = ["sample_id",'file_path'])
df_file_list.to_csv(output_path+'/output/kraken_samples.tsv', sep='\t', header=False, index=False)

