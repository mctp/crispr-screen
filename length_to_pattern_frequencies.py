import gzip
import os
import sys
import argparse
from Bio import SeqIO
import subprocess

parser = argparse.ArgumentParser(description="Filter sequences from a FASTQ file.")
parser.add_argument('--in-fastq', required=True, help="Input FASTQ file (gzipped or plain text).")
parser.add_argument('--pattern', required=True, default="GATCTAGTTACGCC", help="Pattern to search for in sequences.")
parser.add_argument('--nreads', type=int, default=None, help="Number of reads in the FASTQ file to scan (default: all).")
parser.add_argument('--max-tries', type=int, default=10, help="Pattern should be fully captured within this many bases into the sequences (default: 10).")

args = parser.parse_args()
in_fastq = args.in_fastq
nreads = args.nreads
pattern = args.pattern.upper()
max_tries = args.max_tries

if not all(c in 'ACGT' for c in pattern):
    print("Error: Pattern contains invalid characters. Only A, C, G, and T are allowed.")
    sys.exit(1)

def is_available(name):
    try:
        result = subprocess.run([name, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False

def is_gzipped(file_path):
    with open(file_path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def open_fastq(file_path):
    if is_gzipped(file_path):
        if is_available('pigz'):
            return subprocess.Popen(['pigz', '-dc', file_path], stdout=subprocess.PIPE, text=True).stdout
        else:
            return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'rt')

def sliding_pattern_search(fastq, pattern, nreads, max_tries):
    counter = 0
    distance_counts = {}
    representative_sequences = {}
    total_reads = 0
    with open_fastq(fastq) as handle:
        for i, record in enumerate(SeqIO.parse(handle, "fastq")):
            if nreads is not None and i >= nreads:
                break
            total_reads += 1
            sequence = str(record.seq)
            found = False
            for j in range(min(max_tries, len(sequence) - len(pattern) + 1)):
                if sequence[j:j+len(pattern)] == pattern:
                    counter += 1
                    if j in distance_counts:
                        distance_counts[j] += 1
                    else:
                        distance_counts[j] = 1
                    if j not in representative_sequences:
                        representative_sequences[j] = sequence
                    found = True
                    break
            if not found:
                if -1 in distance_counts:
                    distance_counts[-1] += 1
                else:
                    distance_counts[-1] = 1
                if -1 not in representative_sequences:
                    representative_sequences[-1] = sequence

    if total_reads > 0:
        print(f"Pattern found in {counter} ({counter / total_reads * 100:.2f}%) of total reads.")
    sorted_distance_counts = dict(sorted(distance_counts.items(), key=lambda item: item[1], reverse=True))
    print("Distance frequencies and representative sequences:")
    max_length = max_tries + len(pattern) + 25
    max_distance_length = max(len(f"Distance {distance}") for distance in sorted_distance_counts.keys())
    max_count_length = max(len(f"{count}") for count in sorted_distance_counts.values())
    max_percentage_length = max(len(f"{(count / total_reads) * 100:.2f}%") for count in sorted_distance_counts.values())

    for distance, count in sorted_distance_counts.items():
        percentage = (count / total_reads) * 100
        if percentage <= 1:
            continue
        sequence = representative_sequences[distance]
        if len(sequence) > max_length:
            sequence = sequence[:max_length] + "..."
        if distance == -1:
            distance_str = "Distance *"
        else:
            distance_str = f"Distance {distance}"
        
        distance_str = distance_str.ljust(max_distance_length)
        count_str = f"{count}".ljust(max_count_length)
        percentage_str = f"{percentage:.2f}%".ljust(max_percentage_length)
        
        if distance == -1:
            print(f"{distance_str}: {count_str} ({percentage_str}) {sequence:<{max_length + 3}}")
        else:
            pattern_start = sequence.find(pattern)
            if pattern_start != -1:
                sequence = sequence[:pattern_start] + f"\033[91m{pattern}\033[0m" + sequence[pattern_start + len(pattern):]
            print(f"{distance_str}: {count_str} ({percentage_str}) {sequence:<{max_length + 3}}")

sliding_pattern_search(in_fastq, pattern, nreads, max_tries)
