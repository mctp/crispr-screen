# command example: python reorient_fastq.py $US_SEQ $IN_R1_FASTQ $IN_R2_FASTQ $IN_R1_RE_FASTQ $IN_R2_RE_FASTQ $FILE_INTERLEAVED

import os
import subprocess
import shutil
from Bio import SeqIO
import sys
import argparse
import gzip

def parse_args():
    parser = argparse.ArgumentParser(description='Reorient FASTQ files.')
    parser.add_argument('search_sequence', type=str, help='Search sequence')
    parser.add_argument('in_r1_fastq', type=str, help='Input R1 FASTQ file')
    parser.add_argument('in_r2_fastq', type=str, help='Input R2 FASTQ file')
    parser.add_argument('out_r1_re_fastq_gz', type=str, help='Output reoriented R1 FASTQ file')
    parser.add_argument('out_r2_re_fastq_gz', type=str, help='Output reoriented R2 FASTQ file')
    parser.add_argument('file_interleaved_gz', type=str, help='Output interleaved FASTQ file')
    return parser.parse_args()

args = parse_args()

args.search_sequence = args.search_sequence.upper()

# removes .gz from the end of the file name
fastq_interleaved = args.file_interleaved_gz[:-3]
out_r1_re_fastq = args.out_r1_re_fastq_gz[:-3]
out_r2_re_fastq = args.out_r2_re_fastq_gz[:-3]

def interleave(handle_r1, handle_r2, handle_interleaved):
    for fwd, rev in zip(SeqIO.parse(handle_r1, "fastq"), SeqIO.parse(handle_r2, "fastq")):
        assert fwd.id == rev.id
        fwd.id += "/1"
        rev.id += "/2"
        SeqIO.write(fwd, handle_interleaved, "fastq")
        SeqIO.write(rev, handle_interleaved, "fastq")

if not ((os.path.exists(fastq_interleaved) and os.path.getsize(fastq_interleaved) > 0) or 
    (os.path.exists(args.file_interleaved_gz) and os.path.getsize(args.file_interleaved_gz) > 0)):
    with gzip.open(args.in_r1_fastq, "rt") as handle_r1, gzip.open(args.in_r2_fastq, "rt") as handle_r2, open(fastq_interleaved, "w") as handle_interleaved:
        interleave(handle_r1, handle_r2, handle_interleaved)
    print(f"Interleaved FASTQ files written to {fastq_interleaved}")
else:
    print(f"Interleaved FASTQ file {fastq_interleaved} or {args.file_interleaved_gz} already exists and is not empty.")

def reorient_fastq(fastq_interleaved, fastq_r1_out, fastq_r2_out, search_sequence):
    if fastq_interleaved.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
    
    with open_func(fastq_interleaved, mode) as f, open(fastq_r1_out, 'w') as f1, open(fastq_r2_out, 'w') as f2:
        records = SeqIO.parse(f, 'fastq')
        not_written_count = 0
        for fwd in records:
            rev = next(records)
            fwd_min_score = min(fwd.letter_annotations["phred_quality"])
            rev_min_score = min(rev.letter_annotations["phred_quality"])
            # search for the sequence in the forward and reverse reads
            if search_sequence in fwd.seq:
                SeqIO.write(fwd, f1, 'fastq')
                SeqIO.write(rev, f2, 'fastq')
            elif search_sequence in rev.seq:
                SeqIO.write(rev, f1, 'fastq')
                SeqIO.write(fwd, f2, 'fastq')
            else:
                # if the search sequence is not found in either read, write the read with the higher quality score to R1 file
                if fwd_min_score >= 20 and rev_min_score < 20:
                    SeqIO.write(fwd, f1, 'fastq')
                    SeqIO.write(rev, f2, 'fastq')
                elif fwd_min_score < 20 and rev_min_score >= 20:
                    SeqIO.write(rev, f1, 'fastq')
                    SeqIO.write(fwd, f2, 'fastq')
                elif fwd_min_score >= 20 and rev_min_score >= 20:
                    # if both reads have quality scores above 20, write the read with the higher quality score to R1 file
                    if fwd_min_score >= rev_min_score:
                        SeqIO.write(fwd, f1, 'fastq')
                        SeqIO.write(rev, f2, 'fastq')
                    else:
                        SeqIO.write(rev, f1, 'fastq')
                        SeqIO.write(fwd, f2, 'fastq')
                else:
                    # if both reads have quality scores below 20, do not write the reads
                    not_written_count += 1
        print(f"Number of reads not written due to not passing threshold: {not_written_count}")

if  not (os.path.exists(out_r1_re_fastq) and os.path.getsize(out_r1_re_fastq) > 0) or 
    not (os.path.exists(out_r2_re_fastq) and os.path.getsize(out_r2_re_fastq) > 0):
    if os.path.exists(fastq_interleaved) and os.path.getsize(fastq_interleaved) > 0:
        interleaved_file = fastq_interleaved
    elif os.path.exists(args.file_interleaved_gz) and os.path.getsize(args.file_interleaved_gz) > 0:
        interleaved_file = args.file_interleaved_gz
    else:
        print("Error: Interleaved FASTQ file not found.")
        sys.exit(1)
    reorient_fastq(interleaved_file, out_r1_re_fastq, out_r2_re_fastq, args.search_sequence)
    print(f"Reoriented FASTQ files written to {out_r1_re_fastq} and {out_r2_re_fastq}")
else:
    print(f"Reoriented FASTQ files {out_r1_re_fastq} and {out_r2_re_fastq} already exist and are not empty.")


def gzip_file_with_pigz(input_file, output_file, cpus=1):
    try:
        subprocess.run(['pigz', '-p', str(cpus), '-c', input_file], stdout=open(output_file, 'wb'), check=True)
    except FileNotFoundError:
        print("pigz not found. Please install pigz to use this feature.")
        sys.exit(1)

def gzip_file_with_gzip(input_file, output_file, cpus=1):
    # cpus argument is ignored for gzip
    with open(input_file, 'rb') as f_in, gzip.open(output_file, 'wb') as f_out:
        f_out.writelines(f_in)


# Determine which gzip function to use
if shutil.which('pigz'):
    gzip_file = gzip_file_with_pigz
else:
    gzip_file = gzip_file_with_gzip

# def gzip_file(input_file, output_file, cpus=1):
#     try:
#         subprocess.run(['pigz', '-p', str(cpus), '-c', input_file], stdout=open(output_file, 'wb'), check=True)
#     except FileNotFoundError:
#         print("pigz not found. Please install pigz to use this feature.")
#         sys.exit(1)

# def gzip_file(input_file, output_file):
#     with open(input_file, 'rb') as f_in, gzip.open(output_file, 'wb') as f_out:
#         f_out.writelines(f_in)

if os.path.exists(args.out_r1_re_fastq_gz) and os.path.getsize(args.out_r1_re_fastq_gz) > 0:
    print(f"Warning: {args.out_r1_re_fastq_gz} already exists and is not empty.")
else:
    gzip_file(out_r1_re_fastq, args.out_r1_re_fastq_gz,cpus=ncpu)

if os.path.exists(args.out_r2_re_fastq_gz) and os.path.getsize(args.out_r2_re_fastq_gz) > 0:
    print(f"Warning: {args.out_r2_re_fastq_gz} already exists and is not empty.")
else:
    gzip_file(out_r2_re_fastq, args.out_r2_re_fastq_gz,cpus=ncpu)

if os.path.exists(args.file_interleaved_gz) and os.path.getsize(args.file_interleaved_gz) > 0:
    print(f"Warning: {args.file_interleaved_gz} already exists and is not empty.")
else:
    gzip_file(fastq_interleaved, args.file_interleaved_gz,cpus=ncpu)

# Confirm success of gzip process before removing the original unzipped files
gzipped_files = [args.out_r1_re_fastq_gz, args.out_r2_re_fastq_gz, args.file_interleaved_gz]
success = all(os.path.exists(f) and os.path.getsize(f) > 0 for f in gzipped_files)

if success:
    os.remove(out_r1_re_fastq)
    os.remove(out_r2_re_fastq)
    os.remove(fastq_interleaved)
    print(f"Gzipped files written to {args.out_r1_re_fastq_gz}, {args.out_r2_re_fastq_gz}, and {args.file_interleaved_gz}")
else:
    print("Gzip process failed for one or more files. Original files were not removed.")
