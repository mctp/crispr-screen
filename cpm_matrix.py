import argparse
import pandas as pd
import os
from process_metadata import process_metadata

parser = argparse.ArgumentParser(description='Generate count and normalized count matrices.')
parser.add_argument('--config', help='Path to the config file')
parser.add_argument('--metadata', help='Path to the metadata file')
parser.add_argument('--input-fastq-dir', help='Path to the directory containing the FASTQ files')
parser.add_argument('--mode', default="fastq", help='Method used for counting (See config file)')
args = parser.parse_args()

def parse_config(file_path):
    config = {}
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith('#'):
                key_value = line.split('=', 1)
                key = key_value[0].strip()
                value = key_value[1].split('#', 1)[0].strip() if len(key_value) > 1 else ''
                if value.startswith('(') and value.endswith(')'):
                    value = value[1:-1].split()
                config[key] = value
    return config

if args.config is None:
    if not os.path.isfile(args.metadata):
        raise FileNotFoundError(f"Metadata file {args.metadata} does not exist.")
    if not os.path.isdir(args.input_fastq_dir):
        raise NotADirectoryError(f"Input FASTQ directory {args.input_fastq_dir} does not exist or is not a directory.")
    if args.mode is None:
        raise ValueError("Mode is required if config is not provided.")
    os.environ['METADATA_FILE'] = args.metadata
    os.environ['FASTQ_DIR'] = args.input_fastq_dir
    os.environ['MODE'] = args.mode
else:
    config = parse_config(args.config)
    for key, value in config.items():
        if isinstance(value, list):
            value = ' '.join(value)
        os.environ[key] = value

# Expand environment variables in paths
metadata_file = os.path.expandvars(os.environ.get('METADATA_FILE'))
fastq_dir = os.path.expandvars(os.environ.get('FASTQ_DIR'))
output_dir = os.path.expandvars(os.environ.get('OUTPUT_DIR'))

# Debug prints
print(f"METADATA_FILE: {metadata_file}")
print(f"FASTQ_DIR: {fastq_dir}")
print(f"MODE: {os.environ.get('MODE')}")
print(f"OUTPUT_DIR: {output_dir}")

# Print the contents of the metadata file for debugging
with open(metadata_file, 'r') as file:
    print(file.read())

libraries, sample_names, fastq_files = process_metadata(metadata_file, fastq_dir)

count_method = os.environ['MODE']

def compile_sgrna_count_matrix(count_file_paths):
    combined_df = pd.DataFrame()
    for file_path in count_file_paths:
        df = pd.read_csv(file_path, sep='\t', header=0)
        df.rename(columns={'sgRNA': 'sgrna_id', 'Gene': 'sgrna_target'}, inplace=True)
        df.set_index(['sgrna_id', 'sgrna_target'], inplace=True)
        if combined_df.empty:
            combined_df = df
        else:
            combined_df = combined_df.join(df, how='outer')
    combined_df.fillna(0, inplace=True)
    return combined_df

count_file_paths = [os.path.join(output_dir, f'{sample}_{count_method}/{sample}_{count_method}.count.txt') for sample in sample_names]

sgrna_count_matrix = compile_sgrna_count_matrix(count_file_paths)

sgrna_count_matrix.to_csv(f"{output_dir}/sgrna_count_matrix.txt", sep='\t')
sgrna_count_matrix.to_csv(f"{output_dir}/sgrna_count_matrix.csv")

def compute_cpm(count_matrix):
    cpm_matrix = count_matrix.div(count_matrix.sum(axis=0), axis=1) * 1e6
    return cpm_matrix

sgrna_cpm_matrix = compute_cpm(sgrna_count_matrix)
sgrna_cpm_matrix.to_csv(f"{output_dir}/sgrna_cpm_matrix.txt", sep='\t')
sgrna_cpm_matrix.to_csv(f"{output_dir}/sgrna_cpm_matrix.csv")

for sample in sample_names:
    sample_cpm_df = sgrna_cpm_matrix[[sample]].copy()
    sample_cpm_df.columns = ['cpm']
    sample_cpm_df.to_csv(f"{output_dir}/{sample}_{count_method}/{sample}_{count_method}_cpm.txt", sep='\t')
    sample_cpm_df.to_csv(f"{output_dir}/{sample}_{count_method}/{sample}_{count_method}_cpm.csv")
