from Bio import SeqIO
import sys
import argparse
import gzip
import os
from multiprocessing import Pool
import numpy as np
import subprocess

def parse_args():
    parser = argparse.ArgumentParser(description='Reorient FASTQ files.')
    parser.add_argument('--sequences-r1', '--seq-r1', type=str, help='Search sequences for R1 (comma-separated)')
    parser.add_argument('--sequences-r2', '--seq-r2', type=str, help='Search sequences for R2 (comma-separated)')

    parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs to use (default: 1).')
    parser.add_argument('--batch-size', type=int, default=1000, help='Batch size for processing (default: 1000).')
    parser.add_argument('--nreads', type=int, default=None, help='Number of reads to process (default: all).')

    parser.add_argument('--in-fastq-r1', type=str, help='Input R1 FASTQ file (gzipped).')
    parser.add_argument('--in-fastq-r2', type=str, help='Input R2 FASTQ file (gzipped).')
    parser.add_argument('--out-fastq-r1', type=str, help='Output reoriented R1 FASTQ file (gzipped).')
    parser.add_argument('--out-fastq-r2', type=str, help='Output reoriented R2 FASTQ file (gzipped).')
    parser.add_argument('--file-interleaved-gz', type=str, help='Output interleaved FASTQ file (gzipped).')

    return parser.parse_args()

args = parse_args()

ncpu = args.cpus
batch_size = args.batch_size
nreads = args.nreads

search_sequences_r1 = args.sequences_r1.upper().split(',')
search_sequences_r2 = args.sequences_r2.upper().split(',') if args.sequences_r2 else []

if not search_sequences_r1 or search_sequences_r1 == ['']:
    raise ValueError("R1 sequences must not be empty.")
if not search_sequences_r2:
    print("Warning: R2 sequences are empty (allowed).")

print(f"R1 sequences: {search_sequences_r1}")
print(f"R2 sequences: {search_sequences_r2}")

fastq_interleaved = args.file_interleaved_gz[:-3]
out_r1_re_fastq = args.out_fastq_r1[:-3]
out_r2_re_fastq = args.out_fastq_r2[:-3]

pigz_available = False
result = subprocess.run(['pigz', '--version'], capture_output=True, text=True)
if result.returncode == 0:
    print("pigz available.")
    pigz_available = True


def interleave(handle_r1, handle_r2, handle_interleaved, nreads=None):
    count = 0
    chunk_count = 0
    chunk_size = 1000000
    chunk_file = open(f"{handle_interleaved.name}.{chunk_count}", "w")
    
    for fwd, rev in zip(SeqIO.parse(handle_r1, "fastq"), SeqIO.parse(handle_r2, "fastq")):
        if nreads is not None and count >= nreads:
            break
        assert fwd.id == rev.id
        fwd.id += "/1"
        rev.id += "/2"
        SeqIO.write(fwd, handle_interleaved, "fastq")
        SeqIO.write(fwd, chunk_file, "fastq")
        SeqIO.write(rev, handle_interleaved, "fastq")
        SeqIO.write(rev, chunk_file, "fastq")
        count += 1
        
        if count % chunk_size == 0:
            chunk_file.close()
            chunk_count += 1
            chunk_file = open(f"{handle_interleaved.name}.{chunk_count}", "w")
    
    chunk_file.close()

if not ((os.path.exists(fastq_interleaved) and os.path.getsize(fastq_interleaved) > 0) or 
    (os.path.exists(args.file_interleaved_gz) and os.path.getsize(args.file_interleaved_gz) > 0)):
    with gzip.open(args.in_fastq_r1, "rt") as handle_r1, gzip.open(args.in_fastq_r2, "rt") as handle_r2, open(fastq_interleaved, "w") as handle_interleaved:
        interleave(handle_r1, handle_r2, handle_interleaved, nreads=nreads)
    print(f"Interleaved FASTQ files written to {fastq_interleaved}")
else:
    print(f"Interleaved FASTQ file {fastq_interleaved} or {args.file_interleaved_gz} already exists and is not empty.")

# Check if chunked files are present
chunk_count = 0
while os.path.exists(f"{fastq_interleaved}.{chunk_count}"):
    chunk_count += 1

# If chunked files are not all present, generate them by splitting the interleaved file
if chunk_count == 0:
    with open(fastq_interleaved, "r") as handle_interleaved:
        chunk_count = 0
        chunk_size = 1000000
        chunk_file = open(f"{fastq_interleaved}.{chunk_count}", "w")
        count = 0
        for record in SeqIO.parse(handle_interleaved, "fastq"):
            SeqIO.write(record, chunk_file, "fastq")
            count += 1
            if count % chunk_size == 0:
                chunk_file.close()
                chunk_count += 1
                chunk_file = open(f"{fastq_interleaved}.{chunk_count}", "w")
        chunk_file.close()
    print(f"Chunked files generated from {fastq_interleaved}")

def process_batch(batch, search_sequences_r1, search_sequences_r2):
    results = []
    not_written_count = 0
    rejects = []
    count_r1_fwd = 0
    count_r2_rev = 0
    count_r1_rev = 0
    count_r2_fwd = 0

    for fwd, rev in batch:
        # first test for fwd is R1
        if any(seq in fwd.seq for seq in search_sequences_r1):
            results.append((fwd, rev))
            count_r1_fwd += 1
        elif search_sequences_r2 and any(seq in rev.seq for seq in search_sequences_r2):
            results.append((fwd, rev))
            count_r2_rev += 1
        elif any(seq in rev.seq for seq in search_sequences_r1):
            results.append((rev, fwd))
            count_r1_rev += 1
        elif search_sequences_r2 and any(seq in fwd.seq for seq in search_sequences_r2):
            results.append((rev, fwd))
            count_r2_fwd += 1
        else:
            results.append((fwd, rev))
            not_written_count += 1

    stats = {
        'not_written_count': not_written_count,
        'count_r1_fwd': count_r1_fwd,
        'count_r2_rev': count_r2_rev,
        'count_r1_rev': count_r1_rev,
        'count_r2_fwd': count_r2_fwd
    }
    return results, rejects, stats
    

def process_and_write(args):
    batch, search_sequences_r1, search_sequences_r2 = args
    results, rejects, stats = process_batch(batch, search_sequences_r1, search_sequences_r2)
    return results, rejects, stats

def reorient_fastq(fastq_interleaved, fastq_r1_out, fastq_r2_out, search_sequences_r1=search_sequences_r1, search_sequences_r2=search_sequences_r2, batch_size=1000, cpus=1):
    print(f"cpus: {cpus}")
    print(f"batch_size: {batch_size}")
    if fastq_interleaved.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'

    total_stats = {}
    
    with open_func(fastq_interleaved, mode) as f, open(fastq_r1_out, 'w') as f1, open(fastq_r2_out, 'w') as f2, \
         open(fastq_r1_out + '.rejects', 'w') as f1_rejects, open(fastq_r2_out + '.rejects', 'w') as f2_rejects:
        records = SeqIO.parse(f, 'fastq')
        batch = []
        with Pool(processes=cpus) as pool:
            while True:
                try:
                    fwd = next(records)
                    rev = next(records)
                    batch.append((fwd, rev))
                    if len(batch) >= batch_size:
                        batches = [(batch[i::cpus], search_sequences_r1, search_sequences_r2) for i in range(cpus)]
                        results = pool.map(process_and_write, batches)
                        for result in results:
                            res, rejects, stats = result
                            for fwd, rev in res:
                                SeqIO.write(fwd, f1, 'fastq')
                                SeqIO.write(rev, f2, 'fastq')
                            for fwd, rev in rejects:
                                SeqIO.write(fwd, f1_rejects, 'fastq')
                                SeqIO.write(rev, f2_rejects, 'fastq')
                            for key in stats:
                                if key not in total_stats:
                                    total_stats[key] = 0
                                total_stats[key] += stats[key]
                        batch = []
                except StopIteration:
                    if batch:
                        batches = [(batch[i::cpus], search_sequences_r1, search_sequences_r2) for i in range(cpus)]
                        results = pool.map(process_and_write, batches)
                        for result in results:
                            res, not_written_count, rejects, stats = result
                            for fwd, rev in res:
                                SeqIO.write(fwd, f1, 'fastq')
                                SeqIO.write(rev, f2, 'fastq')
                            for fwd, rev in rejects:
                                SeqIO.write(fwd, f1_rejects, 'fastq')
                                SeqIO.write(rev, f2_rejects, 'fastq')
                            for key in stats:
                                if key not in total_stats:
                                    total_stats[key] = 0
                                total_stats[key] += stats[key]
                    break
    # Print summary of total_stats
    for key, value in total_stats.items():
        print(f"{key}: {value}")

if not (os.path.exists(out_r1_re_fastq) and os.path.getsize(out_r1_re_fastq) > 0) or \
   not (os.path.exists(out_r2_re_fastq) and os.path.getsize(out_r2_re_fastq) > 0):
    if os.path.exists(fastq_interleaved) and os.path.getsize(fastq_interleaved) > 0:
        interleaved_file = fastq_interleaved
    elif os.path.exists(args.file_interleaved_gz) and os.path.getsize(args.file_interleaved_gz) > 0:
        interleaved_file = args.file_interleaved_gz
    else:
        print("Error: Interleaved FASTQ file not found.")
        sys.exit(1)
    chunk_count = 0
    while os.path.exists(f"{interleaved_file}.{chunk_count}"):
        chunk_file = f"{interleaved_file}.{chunk_count}"
        reorient_fastq(chunk_file, f"{out_r1_re_fastq}.{chunk_count}", f"{out_r2_re_fastq}.{chunk_count}", search_sequences_r1=search_sequences_r1, search_sequences_r2=search_sequences_r2, cpus=ncpu, batch_size=batch_size)
        chunk_count += 1

    # Combine chunked output files into final output files
    with open(out_r1_re_fastq, 'w') as f1_out, open(out_r2_re_fastq, 'w') as f2_out:
        for i in range(chunk_count):
            with open(f"{out_r1_re_fastq}.{i}", 'r') as f1_in, open(f"{out_r2_re_fastq}.{i}", 'r') as f2_in:
                f1_out.writelines(f1_in.readlines())
                f2_out.writelines(f2_in.readlines())
            os.remove(f"{out_r1_re_fastq}.{i}")
            os.remove(f"{out_r2_re_fastq}.{i}")

    print(f"Reoriented FASTQ files written to {out_r1_re_fastq} and {out_r2_re_fastq}")
else:
    print(f"Reoriented FASTQ files {out_r1_re_fastq} and {out_r2_re_fastq} already exist and are not empty.")

def gzip_file(input_file, output_file):
    if pigz_available:
        subprocess.run(['pigz', '-c', input_file], stdout=open(output_file, 'wb'))
    else:
        with open(input_file, 'rb') as f_in, gzip.open(output_file, 'wb') as f_out:
            f_out.writelines(f_in)

if os.path.exists(args.out_fastq_r1) and os.path.getsize(args.out_fastq_r1) > 0:
    print(f"Warning: {args.out_fastq_r1} already exists and is not empty.")
else:
    gzip_file(out_r1_re_fastq, args.out_fastq_r1)

if os.path.exists(args.out_fastq_r2) and os.path.getsize(args.out_fastq_r2) > 0:
    print(f"Warning: {args.out_fastq_r2} already exists and is not empty.")
else:
    gzip_file(out_r2_re_fastq, args.out_fastq_r2)
