import gzip
import os
import sys
import argparse

parser = argparse.ArgumentParser(description="Filter sequences from a FASTQ file.")
parser.add_argument('--in-fastq', required=True, help="Input FASTQ file (gzipped or plain text).")
parser.add_argument('--pattern', required=True, default="GATCTAGTTACGCC", help="Pattern to search for in sequences.")
parser.add_argument('--max-tries', type=int, default=10, help="Pattern should be fully captured within this many bases into the sequences.")

args = parser.parse_args()
input_file = args.in_fastq
filtered_file = input_file.rsplit('.', 2)[0] + '.seq.gz'
pattern = args.pattern.upper()
max_tries = args.max_tries

if not all(c in 'ACGT' for c in pattern):
    print("Error: Pattern contains invalid characters. Only A, C, G, and T are allowed.")
    sys.exit(1)

def is_gzipped(file_path):
    with open(file_path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def filter_sequences(input_file):
    output_file = input_file.rsplit('.', 2)[0] + '.seq.gz'
    
    if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        print(f"{output_file} already exists and is not empty. Skipping processing.")
        return
    
    if is_gzipped(input_file):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
    
    with open_func(input_file, mode) as fin, gzip.open(output_file, 'wt') as fout:
        for i, line in enumerate(fin):
            if i % 4 == 1: # enumerate starts at 0 hence modulo 4 is 1 to get each second line (sequence)
                fout.write(line)
    
    print(f"Sequences have been written to {output_file}")

def find_pattern_in_sequences(filtered_file, pattern="GATCTAGTTACGCC", max_tries=10):
    counter = 0
    distance_counts = {}

    with gzip.open(filtered_file, 'rt') as f:
        for line in f:
            found = False
            for i in range(max_tries):
                if line[i:i+len(pattern)] == pattern:
                    counter += 1
                    if i in distance_counts:
                        distance_counts[i] += 1
                    else:
                        distance_counts[i] = 1
                    found = True
                    break
            if not found:
                if -1 in distance_counts:
                    distance_counts[-1] += 1
                else:
                    distance_counts[-1] = 1

    print(f"Pattern found {counter} times.")
    print(f"Distance counts: {distance_counts}")


# extract just the sequences from the fastq file
filter_sequences(input_file)

# find pattern and report frequency of distances from start of sequence to pattern
find_pattern_in_sequences(filtered_file, pattern,)
