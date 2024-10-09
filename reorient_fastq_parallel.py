from Bio import SeqIO
import sys
import argparse
import gzip
import os
from multiprocessing import Pool
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Reorient FASTQ files.')
    parser.add_argument('search_sequence', type=str, help='Search sequence')
    parser.add_argument('in_r1_fastq', type=str, help='Input R1 FASTQ file')
    parser.add_argument('in_r2_fastq', type=str, help='Input R2 FASTQ file')
    parser.add_argument('out_r1_re_fastq_gz', type=str, help='Output reoriented R1 FASTQ file')
    parser.add_argument('out_r2_re_fastq_gz', type=str, help='Output reoriented R2 FASTQ file')
    parser.add_argument('file_interleaved_gz', type=str, help='Output interleaved FASTQ file')
    parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs to use')
    parser.add_argument('--batch-size', type=int, default=1000, help='Batch size for processing')
    parser.add_argument('--nreads', type=int, default=None, help='Number of reads to process')
    return parser.parse_args()

args = parse_args()

ncpu = args.cpus
batch_size = args.batch_size
nreads = args.nreads

args.search_sequence = args.search_sequence.upper()

fastq_interleaved = args.file_interleaved_gz[:-3]
out_r1_re_fastq = args.out_r1_re_fastq_gz[:-3]
out_r2_re_fastq = args.out_r2_re_fastq_gz[:-3]

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
    with gzip.open(args.in_r1_fastq, "rt") as handle_r1, gzip.open(args.in_r2_fastq, "rt") as handle_r2, open(fastq_interleaved, "w") as handle_interleaved:
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

def process_batch(batch, search_sequence):
    results = []
    not_written_count = 0
    rejects = []
    for fwd, rev in batch:
        #fwd_median_score = np.median(fwd.letter_annotations["phred_quality"])
        #rev_median_score = np.median(rev.letter_annotations["phred_quality"])
        # ATCTAGTTACGCCAAGC is part of the reverse primer sequence
        if search_sequence in fwd.seq or 'ATCTAGTTACGCCAAGC' in fwd.seq:
            results.append((fwd, rev))
        # 'GGGTAGTTTGCAGTTTT' is part of the forward primer sequence
        elif search_sequence in rev.seq or 'GGGTAGTTTGCAGTTTT' in fwd.seq:
            results.append((rev, fwd))
        else:
            results.append((fwd, rev))
    return results, not_written_count, rejects

def process_and_write(args):
    batch, search_sequence = args
    results, not_written_count, rejects = process_batch(batch, search_sequence)
    return results, not_written_count, rejects

def reorient_fastq(fastq_interleaved, fastq_r1_out, fastq_r2_out, search_sequence, batch_size=1000, cpus=1):
    #debug
    print(f"cpus: {cpus}")
    print(f"batch_size: {batch_size}")
    if fastq_interleaved.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
    
    with open_func(fastq_interleaved, mode) as f, open(fastq_r1_out, 'w') as f1, open(fastq_r2_out, 'w') as f2, \
         open(fastq_r1_out + '.rejects', 'w') as f1_rejects, open(fastq_r2_out + '.rejects', 'w') as f2_rejects:
        records = SeqIO.parse(f, 'fastq')
        batch = []
        total_not_written_count = 0

        with Pool(processes=cpus) as pool:
            while True:
                try:
                    fwd = next(records)
                    rev = next(records)
                    batch.append((fwd, rev))
                    if len(batch) >= batch_size:
                        batches = [(batch[i::cpus], search_sequence) for i in range(cpus)]
                        results = pool.map(process_and_write, batches)
                        for result in results:
                            res, not_written_count, rejects = result
                            for fwd, rev in res:
                                SeqIO.write(fwd, f1, 'fastq')
                                SeqIO.write(rev, f2, 'fastq')
                            for fwd, rev in rejects:
                                SeqIO.write(fwd, f1_rejects, 'fastq')
                                SeqIO.write(rev, f2_rejects, 'fastq')
                            total_not_written_count += not_written_count
                        batch = []
                except StopIteration:
                    if batch:
                        batches = batches = [(batch[i::cpus], search_sequence) for i in range(cpus)]
                        results = pool.map(process_and_write, batches)
                        for result in results:
                            res, not_written_count, rejects = result
                            for fwd, rev in res:
                                SeqIO.write(fwd, f1, 'fastq')
                                SeqIO.write(rev, f2, 'fastq')
                            for fwd, rev in rejects:
                                SeqIO.write(fwd, f1_rejects, 'fastq')
                                SeqIO.write(rev, f2_rejects, 'fastq')
                            total_not_written_count += not_written_count
                    break

        print(f"Number of reads not written due to not passing threshold: {total_not_written_count}")

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
        reorient_fastq(chunk_file, f"{out_r1_re_fastq}.{chunk_count}", f"{out_r2_re_fastq}.{chunk_count}", args.search_sequence, cpus=ncpu, batch_size=batch_size)
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
    with open(input_file, 'rb') as f_in, gzip.open(output_file, 'wb') as f_out:
        f_out.writelines(f_in)

if os.path.exists(args.out_r1_re_fastq_gz) and os.path.getsize(args.out_r1_re_fastq_gz) > 0:
    print(f"Warning: {args.out_r1_re_fastq_gz} already exists and is not empty.")
else:
    gzip_file(out_r1_re_fastq, args.out_r1_re_fastq_gz)

if os.path.exists(args.out_r2_re_fastq_gz) and os.path.getsize(args.out_r2_re_fastq_gz) > 0:
    print(f"Warning: {args.out_r2_re_fastq_gz} already exists and is not empty.")
else:
    gzip_file(out_r2_re_fastq, args.out_r2_re_fastq_gz)

# if os.path.exists(args.file_interleaved_gz) and os.path.getsize(args.file_interleaved_gz) > 0:
#     print(f"Warning: {args.file_interleaved_gz} already exists and is not empty.")
# else:
#     gzip_file(fastq_interleaved, args.file_interleaved_gz)

# gzipped_files = [args.out_r1_re_fastq_gz, args.out_r2_re_fastq_gz, args.file_interleaved_gz]
# success = all(os.path.exists(f) and os.path.getsize(f) > 0 for f in gzipped_files)

# if success:
#     os.remove(out_r1_re_fastq)
#     os.remove(out_r2_re_fastq)
#     os.remove(fastq_interleaved)
#     print(f"Gzipped files written to {args.out_r1_re_fastq_gz}, {args.out_r2_re_fastq_gz}, and {args.file_interleaved_gz}")
# else:
#     print("Gzip process failed for one or more files. Original files were not removed.")
