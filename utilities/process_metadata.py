import os
import yaml
import pandas as pd

def process_metadata(metadata_file, fastq_dir=""):
    libraries = []
    if (metadata_file.endswith('.yml') or metadata_file.endswith('.yaml')):
        with open(metadata_file, 'r') as f:
            data = yaml.safe_load(f)
        
        # check the structure
        if 'samples' not in data:
            raise ValueError(f"YAML file must contain 'samples' key. Found keys: {list(data.keys())}")
        samples = data['samples']
        if not isinstance(samples, list):
            raise ValueError("'samples' must be a list in YAML file")
        
        libraries = []
        sample_names = []
        fastq_files = []
        
        for sample in samples:
               
            # Validate required fields
            if not isinstance(sample, dict):
                raise ValueError(f"Each sample must be a dictionary, got: {type(sample)}")

            # Skip commented samples
            if sample.get('sample', '').startswith('#'):
                continue

            if 'library' not in sample or 'sample' not in sample:
                raise ValueError(f"Each sample must have 'library' and 'sample' fields. Got: {list(sample.keys())}")
            
            libraries.append(sample['library'])
            sample_names.append(sample['sample'])
            
            r1_files = sample.get('fastq_r1', [])
            r2_files = sample.get('fastq_r2', [])
            
            # Ensure they are lists
            if isinstance(r1_files, str):
                r1_files = [r1_files]
            if isinstance(r2_files, str):
                r2_files = [r2_files]
            
            # Add directory prefix if specified
            if fastq_dir:
                r1_files = [os.path.join(fastq_dir, f) for f in r1_files]
                r2_files = [os.path.join(fastq_dir, f) for f in r2_files]
            
            # Format as comma-separated with semicolon separator
            if r2_files:
                fastq_files.append(",".join(r1_files) + ";" + ",".join(r2_files))
            else:
                fastq_files.append(",".join(r1_files))
                
    else:
        # handle sample metadata in txt format
        df = pd.read_csv(metadata_file, delimiter='\t', comment='#')
        libraries = df.iloc[:, 0].tolist()
        sample_names = df.iloc[:, 1].tolist()
        fastq_files = []

        for _, row in df.iterrows():
            fq1s = row.iloc[2].split(',')
            fq2s = row.iloc[3].split(',') if pd.notna(row.iloc[3]) and row.iloc[3] else []
            fq1s = [os.path.join(fastq_dir, fq1) for fq1 in fq1s]
            fq2s = [os.path.join(fastq_dir, fq2) for fq2 in fq2s]
            
            if fq2s:
                fastq_files.append(",".join(fq1s) + ";" + ",".join(fq2s))
            else:
                fastq_files.append(",".join(fq1s))

    return libraries, sample_names, fastq_files
