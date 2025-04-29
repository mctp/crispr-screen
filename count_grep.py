import os
import gzip
import subprocess
from multiprocessing import Pool
import argparse

def revcomp(seq):
    complement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def unzip_and_extract_sequences(fastq_file, output_file):
    with gzip.open(fastq_file, 'rt') as f_in, open(output_file, 'w') as f_out:
        for i, line in enumerate(f_in):
            if i % 4 == 1:
                f_out.write(line)

def count_sgrna(args):
    sgrna_id, sgrna_seq, sgrna_target, in_r1_fastq_seq, in_r2_fastq_seq, output_file = args
    r1_count = subprocess.check_output(['grep', '-c', sgrna_seq, in_r1_fastq_seq]).decode().strip()
    r2_count = subprocess.check_output(['grep', '-c', sgrna_seq, in_r2_fastq_seq]).decode().strip()
    r1_rc_count = subprocess.check_output(['grep', '-c', revcomp(sgrna_seq), in_r1_fastq_seq]).decode().strip()
    r2_rc_count = subprocess.check_output(['grep', '-c', revcomp(sgrna_seq), in_r2_fastq_seq]).decode().strip()
    with open(output_file, 'a') as f:
        f.write(f"{sgrna_id}\t{sgrna_seq}\t{sgrna_target}\t{r1_count}\t{r2_count}\t{r1_rc_count}\t{r2_rc_count}\n")

def count_grep(sample, mode, output_dir, in_r1_fastq, in_r2_fastq, sgrna_list_name, ncpu):
    if in_r1_fastq.endswith('.fq.gz'):
        in_r1_fastq_seq = f"{in_r1_fastq[:-6]}.txt"
    elif in_r1_fastq.endswith('.fastq.gz'):
        in_r1_fastq_seq = f"{in_r1_fastq[:-9]}.txt"
    else:
        raise ValueError("Unsupported file extension for in_r1_fastq")

    if in_r2_fastq.endswith('.fq.gz'):
        in_r2_fastq_seq = f"{in_r2_fastq[:-6]}.txt"
    elif in_r2_fastq.endswith('.fastq.gz'):
        in_r2_fastq_seq = f"{in_r2_fastq[:-9]}.txt"
    else:
        raise ValueError("Unsupported file extension for in_r2_fastq")
    output_file = f"{output_dir}/{sample}_{mode}/{sample}_count.txt"

    if os.path.exists(output_file):
        os.remove(output_file)

    if not os.path.exists(in_r1_fastq_seq):
        print("R1: unzipping fastq and keeping sequences...")
        unzip_and_extract_sequences(in_r1_fastq, in_r1_fastq_seq)

    if not os.path.exists(in_r2_fastq_seq):
        print("R2: unzipping fastq and keeping sequences...")
        unzip_and_extract_sequences(in_r2_fastq, in_r2_fastq_seq)

    print(f"sample: {sample}")
    sgrna_list_file = f"{sgrna_list_name}.txt"
    print(f"counting sgRNAs in {sgrna_list_file}...")

    with open(sgrna_list_file) as f:
        lines = f.readlines()[1:]

    sgrna_data = [(line.split('\t')[0], line.split('\t')[1], line.split('\t')[2].strip(), in_r1_fastq_seq, in_r2_fastq_seq, output_file) for line in lines]

    os.makedirs(f"{output_dir}/{sample}_{mode}", exist_ok=True)
    with open(output_file, 'w') as f:
        f.write("sgrna_id\tsgrna_sequence\tsgrna_target\tR1_count\tR2_count\tR1_rc_count\tR2_rc_count\n")

    with Pool(ncpu) as pool:
        pool.map(count_sgrna, sgrna_data)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count sgRNAs in FASTQ files.')
    parser.add_argument('--sample', required=True, help='Sample name')
    parser.add_argument('--mode', required=True, help='Mode of operation')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    parser.add_argument('--in_r1_fastq', required=True, help='Input R1 FASTQ file')
    parser.add_argument('--in_r2_fastq', required=True, help='Input R2 FASTQ file')
    parser.add_argument('--sgrna_list_name', required=True, help='sgRNA list name')
    parser.add_argument('--ncpu', type=int, default=1, help='Number of CPUs to use')
    args = parser.parse_args()
    count_grep(
        sample=args.sample,
        mode=args.mode,
        output_dir=args.output_dir,
        in_r1_fastq=args.in_r1_fastq,
        in_r2_fastq=args.in_r2_fastq,
        sgrna_list_name=args.sgrna_list_name,
        ncpu=args.ncpu
    )
