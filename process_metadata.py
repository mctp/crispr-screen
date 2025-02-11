import os
import sys
import yaml
import pandas as pd

def process_metadata(metadata_file, fastq_dir=""):
    libraries = []  # Initialize libraries to an empty list
    if (metadata_file.endswith('.yml') or metadata_file.endswith('.yaml')):
        with open(metadata_file, 'r') as f:
            data = yaml.safe_load(f)
        libraries = [sample['library'] for sample in data['samples'] if not sample['sample'].startswith('#')]
        sample_names = [sample['sample'] for sample in data['samples'] if not sample['sample'].startswith('#')]
        fastq_files = [
            ",".join(sample['fastq_r1']) + ";" + ",".join(sample['fastq_r2'])
            for sample in data['samples'] if not sample['sample'].startswith('#')
        ]
    else:
        df = pd.read_csv(metadata_file, delimiter='\t', comment='#')  # Specify the delimiter and ignore commented lines
        libraries = df.iloc[:, 0].tolist()
        sample_names = df.iloc[:, 1].tolist()
        fastq_files = []

        for _, row in df.iterrows():
            fq1s = row.iloc[2].split(',')
            fq2s = row.iloc[3].split(',')
            fq1s = [os.path.join(fastq_dir, fq1) for fq1 in fq1s]
            fq2s = [os.path.join(fastq_dir, fq2) for fq2 in fq2s]
            fastq_files.append(",".join(fq1s) + ";" + ",".join(fq2s))

    return libraries, sample_names, fastq_files

# Example usage:
# libraries, sample_names, fastq_files = process_metadata('path/to/metadata_file', 'path/to/fastq_dir')
