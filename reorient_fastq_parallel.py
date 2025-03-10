from Bio import SeqIO
import sys
import argparse
import gzip
import os
import re
from multiprocessing import Pool
import subprocess
import read_reorientation_barplot as plt
import logging
import parse_config as pc
import shutil
import time


logging.basicConfig(
    level=logging.INFO,
    format="[{asctime}] [{levelname}] {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("output/reorient_fastq_parallel.log")
    ]
)

# check if docker environment
# note: this method may not always be reliable
def is_docker_env():
    try:
        with open('/proc/1/cgroup', 'rt') as f:
            return 'docker' in f.read()
    except FileNotFoundError:
        return False

def check_pigz():
    try:
        result = subprocess.run(["pigz", "--version"], capture_output=True, text=True)
        if result.returncode == 0:
            version_line = result.stdout.split('\n')[0]
            version = version_line.split()[1]
            logging.info(f"pigz version: {version}")
            return True
        else:
            raise FileNotFoundError
    except FileNotFoundError:
        logging.warning("pigz is not available. Falling back to gzip.")
        return False

def check_fastq_file(file_path):
    if not os.path.exists(file_path):
        logging.error(f"File {file_path} does not exist.")
        raise FileNotFoundError(f"File {file_path} does not exist.")
    if os.path.getsize(file_path) == 0:
        logging.error(f"File {file_path} is empty.")
        raise ValueError(f"File {file_path} is empty.")
    try:
        with gzip.open(file_path, "rt") as handle:
            SeqIO.parse(handle, "fastq")
    except Exception as e:
        logging.error(f"File {file_path} is not a valid FASTQ file: {e}")
        raise ValueError(f"File {file_path} is not a valid FASTQ file: {e}")

def interleave(handle_r1, handle_r2, handle_interleaved, nreads=None):
    count = 0
    for record_r1, record_r2 in zip(SeqIO.parse(handle_r1, "fastq"), SeqIO.parse(handle_r2, "fastq")):
        SeqIO.write(record_r1, handle_interleaved, "fastq")
        SeqIO.write(record_r2, handle_interleaved, "fastq")
        count += 1
        if nreads and count >= nreads:
            break

def interleave_seqtk(in_r1, in_r2, out_interleaved_file, nreads=None):
    try:
        if nreads:
            cmd = f"seqtk mergepe <(seqtk seq {in_r1} | head -n {nreads * 4}) <(seqtk seq {in_r2} | head -n {nreads * 4})"
            with open(out_interleaved_file, 'w') as f:
                subprocess.run(cmd, shell=True, executable='/bin/bash', stdout=f, check=True)
        else:
            with open(out_interleaved_file, 'w') as f:
                subprocess.run(['seqtk', 'mergepe', in_r1, in_r2], stdout=f, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"seqtk mergepe failed: {e}")
        raise

def split_interleaved_file(fastq_interleaved, chunk_size=1000000):
    with open(fastq_interleaved, "r") as handle_interleaved:
        chunk_count = 0
        logging.info(f"Writing chunk file {chunk_count} ...")
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
    return chunk_count

def split_interleaved_file_seqtk(fastq_interleaved, chunk_size=1000000):
    try:
        chunk_count = 0
        with open(fastq_interleaved, "r") as f:
            for i, _ in enumerate(f):
                pass
        total_chunks = (i + 1 + chunk_size * 8 - 1) // (chunk_size * 8)
        for chunk in range(total_chunks):
            start_line = chunk * chunk_size * 8
            chunk_file = f"{fastq_interleaved}.{chunk}"
            cmd = f"seqtk seq {fastq_interleaved} | tail -n +{start_line + 1} | head -n {chunk_size * 8} > {chunk_file}"
            logging.info(f"Writing chunk file {chunk} ...")
            # logging.info(f"Running command: {cmd}")
            subprocess.run(cmd, shell=True, executable='/bin/bash', check=True)
            chunk_count += 1
    except subprocess.CalledProcessError as e:
        logging.error(f"seqtk seq failed: {e}")
        raise
    return chunk_count

def process_batch(batch, search_sequences_r1, search_sequences_r2):
    results = []
    count_not_found = 0
    count_r1_fwd = 0
    count_r2_rev = 0
    count_r1_rev = 0
    count_r2_fwd = 0

    for fwd, rev in batch:
        fwd.id = fwd.id + "/1"
        rev.id = rev.id + "/2"
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
            count_not_found += 1

    stats = {
        'count_not_found': count_not_found,
        'count_r1_fwd': count_r1_fwd,
        'count_r2_rev': count_r2_rev,
        'count_r1_rev': count_r1_rev,
        'count_r2_fwd': count_r2_fwd
    }
    return results, stats

def gzip_file(input_file):
    temp_file = f"{input_file}.tmp"
    shutil.copy(input_file, temp_file)
    if pigz_available:
        gzip = "pigz"
    else:
        gzip = "gzip"
    gzip_command = f"{gzip} {temp_file}"
    command = f"nohup sh -c '{gzip_command} && mv {temp_file}.gz {input_file}.gz' >/dev/null 2>&1 &"
    process = subprocess.Popen(command, shell=True)
    return process

def process_and_write(args):
    batch, search_sequences_r1, search_sequences_r2 = args
    results, stats = process_batch(batch, search_sequences_r1, search_sequences_r2)
    return results, stats

def get_fastq_ids(file_paths):
    ids = set()
    for file_path in file_paths:
        if os.path.isfile(file_path):
            with gzip.open(file_path, "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    ids.add(record.id)
        else:
            logging.error(f"Error: {file_path} is not a valid file.")
            raise ValueError(f"{file_path} is not a valid file.")
    return ids

def strip_suffix(read_id):
    return re.sub(r'/[12]$', '', read_id)

def validate_fastq(fastq_files_orig, fastq_files_re):

    # note that this function wants a list of file paths, not a single file path, so we need to wrap the single file paths in brackets
    ids_original_r1 = get_fastq_ids([fastq_files_orig[0]])
    ids_original_r2 = get_fastq_ids([fastq_files_orig[1]])
    ids_reoriented_r1 = get_fastq_ids([fastq_files_re[0]])
    ids_reoriented_r2 = get_fastq_ids([fastq_files_re[1]])

    logging.info(f"Total read IDs in {fastq_files_orig[0]}: {len(ids_original_r1)}")
    logging.info(f"Total read IDs in {fastq_files_orig[1]}: {len(ids_original_r2)}")
    logging.info(f"Total read IDs in {fastq_files_re[0]}: {len(ids_reoriented_r1)}")
    logging.info(f"Total read IDs in {fastq_files_re[1]}: {len(ids_reoriented_r2)}")

    # strip /1 /2 etc to compare read IDs
    ids_original_r1_stripped = {strip_suffix(id) for id in ids_original_r1}
    ids_original_r2_stripped = {strip_suffix(id) for id in ids_original_r2}
    ids_reoriented_r1_stripped = {strip_suffix(id) for id in ids_reoriented_r1}
    ids_reoriented_r2_stripped = {strip_suffix(id) for id in ids_reoriented_r2}

    # compare stripped read IDs
    if ids_original_r1_stripped != ids_reoriented_r1_stripped or ids_original_r2_stripped != ids_reoriented_r2_stripped:
        logging.error("FASTQ files do not contain the same read IDs.")
        raise ValueError("FASTQ files do not contain the same read IDs.")

    logging.info("Read IDs in original and reoriented FASTQ files match.")

    count_reoriented_r1_end1 = sum(1 for id in ids_reoriented_r1 if id.endswith('/1'))
    count_reoriented_r1_end2 = sum(1 for id in ids_reoriented_r1 if id.endswith('/2'))
    count_reoriented_r2_end1 = sum(1 for id in ids_reoriented_r2 if id.endswith('/1'))
    count_reoriented_r2_end2 = sum(1 for id in ids_reoriented_r2 if id.endswith('/2'))

    logging.info(f"Total /1 read IDs in reoriented R1: {count_reoriented_r1_end1}, R2: {count_reoriented_r2_end1}")
    logging.info(f"Total /2 read IDs in reoriented R1: {count_reoriented_r1_end2}, R2: {count_reoriented_r2_end2}") 

def reorient_fastq(fastq_interleaved, fastq_r1_out, fastq_r2_out, search_sequences_r1, search_sequences_r2, batch_size=1000, cpus=1):
    if fastq_interleaved.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'

    total_stats = {}
    
    with open_func(fastq_interleaved, mode) as f, open(fastq_r1_out, 'w') as f1, open(fastq_r2_out, 'w') as f2:
        records = SeqIO.parse(f, 'fastq')
        batch = []
        batch_count = 0
        with Pool(processes=cpus) as pool:
            while True:
                try:
                    fwd = next(records)
                    rev = next(records)
                    batch.append((fwd, rev))
                    if len(batch) >= batch_size:
                        batch_count += 1
                        batches = [(batch[i::cpus], search_sequences_r1, search_sequences_r2) for i in range(cpus)]
                        results = pool.map(process_and_write, batches)
                        for result in results:
                            res, stats = result
                            for fwd, rev in res:
                                SeqIO.write(fwd, f1, 'fastq')
                                SeqIO.write(rev, f2, 'fastq')
                            for key in stats:
                                if key not in total_stats:
                                    total_stats[key] = 0
                                total_stats[key] += stats[key]
                        batch = []
                except StopIteration:
                    if batch:
                        batch_count += 1
                        logging.info(f"Processing final batch {batch_count} ...")
                        batches = [(batch[i::cpus], search_sequences_r1, search_sequences_r2) for i in range(cpus)]
                        results = pool.map(process_and_write, batches)
                        for result in results:
                            res, stats = result
                            for fwd, rev in res:
                                SeqIO.write(fwd, f1, 'fastq')
                                SeqIO.write(rev, f2, 'fastq')
                            for key in stats:
                                if key not in total_stats:
                                    total_stats[key] = 0
                                total_stats[key] += stats[key]
                    break
    return total_stats

def fibonacci_sequence(n):
    fib_seq = [1, 1]
    for _ in range(2, n):
        fib_seq.append(fib_seq[-1] + fib_seq[-2])
    return fib_seq

def check_file_exists(file_path, description):
    fib_seq = fibonacci_sequence(10)  # Pre-fetch 10 Fibonacci values
    for sleep_time in fib_seq:
        if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
            return
        if sleep_time > 10:
            logging.warning(f"Warning: {description} {file_path} does not exist or is 0 bytes. Retrying in {sleep_time} seconds...")
        time.sleep(sleep_time)
    logging.error(f"Error: {description} {file_path} does not exist or is 0 bytes after multiple attempts.")
    sys.exit(1)
def main():
    global pigz_available
    pigz_available = check_pigz()
    global is_docker
    is_docker = is_docker_env()

    parser = argparse.ArgumentParser(description='Reorient FASTQ files.')
    parser.add_argument('--config', type=str, help='Configuration file.')

    parser.add_argument('--sequences-r1', '--seq-r1', type=str, required=True, help='Search sequences for R1 (comma-separated)')
    parser.add_argument('--sequences-r2', '--seq-r2', type=str, required=True, help='Search sequences for R2 (comma-separated)')

    parser.add_argument('--cpus','-p', type=int, default=1, help='Number of CPUs to use (default: 1).')
    parser.add_argument('--batch-size', '-b', type=int, default=1000, help='Batch size for processing (default: 1000).')
    parser.add_argument('--chunk-size', '-c', type=int, default=1000000, help='Chunk size for splitting interleaved file (default: 1000000).')
    parser.add_argument('--nreads', '-n', type=int, default=None, help='Number of reads to process (default: all).')

    parser.add_argument('--in-fastq-r1', '--in-r1', type=str, required=True, help='Input R1 FASTQ file (gzipped).')
    parser.add_argument('--in-fastq-r2', '--in-r2', type=str, required=True, help='Input R2 FASTQ file (gzipped).')
    parser.add_argument('--out-fastq-r1', '--out-r1', type=str, required=True, help='Output reoriented R1 FASTQ file (gzipped).')
    parser.add_argument('--out-fastq-r2', '--out-r2', type=str, required=True, help='Output reoriented R2 FASTQ file (gzipped).')

    parser.add_argument('--file-interleaved-gz', '--int', '--interleaved-fastq', type=str, help='Output interleaved FASTQ file (gzipped).')

    parser.add_argument('--plot', action='store_true', help='Plot sequence type proportions.')
    parser.add_argument('--plot-prefix', type=str, help='Prefix for plot files.')
    parser.add_argument('--plot-search-sequence', '--plot-search-sequences', '--plot-seq', type=str, help='Search sequence(s) for plotting.')
    parser.add_argument('--plot-only', action='store_true', help='Skip read reorientation and do plots.')

    parser.add_argument('--validate', action='store_true', help='Validate FASTQ files after processing.')

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if not args.in_fastq_r1 or not args.in_fastq_r2 or not args.out_fastq_r1 or not args.out_fastq_r2:
        logging.error("Error: Missing required FASTQ file arguments.")
        sys.exit(1)
        config = pc.parse_config(args.config)
        if is_docker:
            output_dir = config['DOCKER_OUTPUT_DIR']
        else:
            output_dir = config['OUTPUT_DIR']
        
    search_sequences_r1 = args.sequences_r1.upper().split(',')
    search_sequences_r2 = args.sequences_r2.upper().split(',') if args.sequences_r2 else []

    if not args.plot and args.plot_only:
        logging.error("Error: --plot-only specified without --plot.")
        sys.exit(1)

    if not args.plot_only:

        logging.info(f"Checking input FASTQ file {args.in_fastq_r1}")
        check_fastq_file(args.in_fastq_r1)
        logging.info(f"Checking input FASTQ file {args.in_fastq_r2}")
        check_fastq_file(args.in_fastq_r2)

        ncpu = args.cpus
        batch_size = args.batch_size
        nreads = args.nreads

        if not search_sequences_r1 or search_sequences_r1 == ['']:
            raise ValueError("R1 sequences must not be empty.")
        if not search_sequences_r2:
            logging.warning("R2 sequences are empty (allowed).")

        logging.info(f"R1 search sequences: {search_sequences_r1}")
        logging.info(f"R2 search sequences: {search_sequences_r2}")

        if args.file_interleaved_gz:
            fastq_interleaved = re.sub(r'\.gz$', '', args.file_interleaved_gz)
        else:
            if "R1" in args.in_fastq_r1:
                fastq_interleaved = re.sub(r'R1', 'interleaved', args.in_fastq_r1)
                fastq_interleaved = re.sub(r'\.(fastq|fq)\.gz$', '.fq', fastq_interleaved)
            else:
                r1_prefix = re.sub(r'\.(fastq|fq)(\.gz)?$', '', args.in_fastq_r1)
                fastq_interleaved = f"{r1_prefix}_interleaved.fq"
            args.file_interleaved_gz = f"{fastq_interleaved}.gz"

        # Verify that args.file_interleaved_gz is identical to fastq_interleaved with ".gz" appended
        if args.file_interleaved_gz != f"{fastq_interleaved}.gz":
            logging.error("Error: args.file_interleaved_gz must be identical to fastq_interleaved with '.gz' appended.")
            sys.exit(1)

        out_r1_re_fastq = args.out_fastq_r1[:-3]
        out_r2_re_fastq = args.out_fastq_r2[:-3]

        # Check if output files already exist and are not empty
        if os.path.exists(args.out_fastq_r1) and os.path.getsize(args.out_fastq_r1) > 0 and \
            os.path.exists(args.out_fastq_r2) and os.path.getsize(args.out_fastq_r2) > 0:
            logging.info(f"Final output files {args.out_fastq_r1} and {args.out_fastq_r2} already exist and are not empty. Skipping processing.")
            
        else:
            if not ((os.path.exists(fastq_interleaved) and os.path.getsize(fastq_interleaved) > 0) or \
                (os.path.exists(args.file_interleaved_gz) and os.path.getsize(args.file_interleaved_gz) > 0)):
                # with gzip.open(args.in_fastq_r1, "rt") as handle_r1, gzip.open(args.in_fastq_r2, "rt") as handle_r2, open(fastq_interleaved, "w") as handle_interleaved:
                logging.info(f"Interleaving FASTQ files {args.in_fastq_r1} and {args.in_fastq_r2} (using seqtk) ...")
                interleave_seqtk(args.in_fastq_r1, args.in_fastq_r2, fastq_interleaved, nreads=nreads)
                logging.info(f"Interleaved FASTQ written to {fastq_interleaved}")
            else:
                logging.info(f"Interleaved FASTQ file {fastq_interleaved} or {args.file_interleaved_gz} already exists and is not empty.")

            # gzip interleaved file if the gz version does not exist
            if not os.path.exists(args.file_interleaved_gz) or os.path.getsize(args.file_interleaved_gz) == 0:
                logging.info(f"Compressing interleaved FASTQ file to {args.file_interleaved_gz}")
                gzip_interleaved_process = gzip_file(fastq_interleaved)

            # Check if unzipped interleaved file exists and if not, unzip
            if not os.path.exists(fastq_interleaved) or os.path.getsize(fastq_interleaved) == 0:
                if os.path.exists(args.file_interleaved_gz) and os.path.getsize(args.file_interleaved_gz) > 0:
                    logging.info(f"Unzipping {args.file_interleaved_gz} to {fastq_interleaved}")
                    with gzip.open(args.file_interleaved_gz, 'rt') as f_in, open(fastq_interleaved, 'w') as f_out:
                        f_out.writelines(f_in)
                else:
                    logging.error(f"Error: Interleaved gzipped file {args.file_interleaved_gz} does not exist or is empty.")
                    sys.exit(1)

            # Check if chunked files are present
            chunk_count = 0
            while os.path.exists(f"{fastq_interleaved}.{chunk_count}"):
                chunk_count += 1

            # If chunked files are not all present, generate them by splitting the interleaved file
            if chunk_count == 0:
                chunk_count = split_interleaved_file_seqtk(fastq_interleaved, chunk_size=args.chunk_size)
            logging.info(f"{chunk_count} chunked files generated from {fastq_interleaved}")

            if not (os.path.exists(out_r1_re_fastq) and os.path.getsize(out_r1_re_fastq) > 0) or \
            not (os.path.exists(out_r2_re_fastq) and os.path.getsize(out_r2_re_fastq) > 0):
                if os.path.exists(fastq_interleaved) and os.path.getsize(fastq_interleaved) > 0:
                    interleaved_file = fastq_interleaved
                elif os.path.exists(args.file_interleaved_gz) and os.path.getsize(args.file_interleaved_gz) > 0:
                    interleaved_file = args.file_interleaved_gz
                else:
                    logging.error("Error: Interleaved FASTQ file not found.")
                    sys.exit(1)
                chunk_count = 0
                logging.info(f"cpus: {ncpu}")
                logging.info(f"batch size: {batch_size}")
                total_chunk_count = 0
                while os.path.exists(f"{interleaved_file}.{chunk_count}"):
                    chunk_count += 1
                    total_chunk_count += 1
                chunk_count = 0
                while os.path.exists(f"{interleaved_file}.{chunk_count}"):
                    logging.info(f"Reorienting chunk file (.{chunk_count}) {chunk_count + 1 } of {total_chunk_count}...")
                    chunk_file = f"{interleaved_file}.{chunk_count}"
                    out_r1_chunk = f"{out_r1_re_fastq}.{chunk_count}"
                    out_r2_chunk = f"{out_r2_re_fastq}.{chunk_count}"

                    if not (os.path.exists(out_r1_chunk) and os.path.getsize(out_r1_chunk) > 0) or \
                    not (os.path.exists(out_r2_chunk) and os.path.getsize(out_r2_chunk) > 0):
                        total_stats = reorient_fastq(chunk_file, out_r1_chunk, out_r2_chunk, search_sequences_r1=search_sequences_r1, search_sequences_r2=search_sequences_r2, cpus=ncpu, batch_size=batch_size)
                        # Print summary of total_stats
                        for key, value in total_stats.items():
                            logging.info(f"{key}: {value}")
                    else:
                        logging.info(f"Reoriented chunk files {out_r1_chunk} and {out_r2_chunk} already exist and are not empty.")

                    chunk_count += 1
                # Combine chunked output files into final output files
                with open(out_r1_re_fastq, 'w') as f1_out, open(out_r2_re_fastq, 'w') as f2_out:
                    for i in range(chunk_count):
                        with open(f"{out_r1_re_fastq}.{i}", 'r') as f1_in, open(f"{out_r2_re_fastq}.{i}", 'r') as f2_in:
                            f1_out.writelines(f1_in.readlines())
                            f2_out.writelines(f2_in.readlines())
                        # os.remove(f"{out_r1_re_fastq}.{i}")
                        # os.remove(f"{out_r2_re_fastq}.{i}")

                logging.info(f"Reoriented FASTQ files written to {out_r1_re_fastq} and {out_r2_re_fastq}")
            else:
                logging.info(f"Reoriented FASTQ files {out_r1_re_fastq} and {out_r2_re_fastq} already exist and are not empty.")

            if os.path.exists(args.out_fastq_r1) and os.path.getsize(args.out_fastq_r1) > 0:
                logging.warning(f"Warning: {args.out_fastq_r1} already exists and is not empty.")
            else:
                logging.info(f"Compressing {out_r1_re_fastq} to {args.out_fastq_r1}")
                gzip_out_r1_process = gzip_file(out_r1_re_fastq)

            if os.path.exists(args.out_fastq_r2) and os.path.getsize(args.out_fastq_r2) > 0:
                logging.warning(f"Warning: {args.out_fastq_r2} already exists and is not empty.")
            else:
                logging.info(f"Compressing {out_r2_re_fastq} to {args.out_fastq_r2}")
                gzip_out_r2_process = gzip_file(out_r2_re_fastq)

            # # Remove auto-incremented chunk files and any non-gzipped FASTQ files
            # for i in range(chunk_count):
            #     chunk_file = f"{fastq_interleaved}.{i}"
            #     if os.path.exists(chunk_file):
            #         os.remove(chunk_file)
            #     out_r1_chunk = f"{out_r1_re_fastq}.{i}"
            #     if os.path.exists(out_r1_chunk):
            #         os.remove(out_r1_chunk)
            #     out_r2_chunk = f"{out_r2_re_fastq}.{i}"
            #     if os.path.exists(out_r2_chunk):
            #         os.remove(out_r2_chunk)

            # if os.path.exists(fastq_interleaved):
            #     os.remove(fastq_interleaved)
            # if os.path.exists(out_r1_re_fastq):
            #     os.remove(out_r1_re_fastq)
            # if os.path.exists(out_r2_re_fastq):
            #     os.remove(out_r2_re_fastq)
    else:
        logging.info("Skipping read reorientation.")
        # Check all 4 input FASTQ files for existence and size > 0 bytes
        input_fastq_files = [args.in_fastq_r1, args.in_fastq_r2, args.out_fastq_r1, args.out_fastq_r2]
        for fastq_file in input_fastq_files:
            if not os.path.exists(fastq_file):
                logging.error(f"Error: Input FASTQ file {fastq_file} does not exist.")
                sys.exit(1)
            if os.path.getsize(fastq_file) == 0:
                logging.error(f"Error: Input FASTQ file {fastq_file} is empty.")
                sys.exit(1)

    if 'gzip_out_r1_process' in locals() or 'gzip_out_r2_process' in locals() or 'gzip_interleaved_process' in locals():
        logging.info("Waiting for gzip processes to finish...")
        if 'gzip_out_r1_process' in locals():
            if gzip_out_r1_process.wait() != 0:
                logging.error("gzip process for R1 failed.")
                sys.exit(1)
            logging.info("    R1 finished...")
        if 'gzip_out_r2_process' in locals():
            if gzip_out_r2_process.wait() != 0:
                logging.error("gzip process for R2 failed.")
                sys.exit(1)
            logging.info("    R2 finished...")
        if 'gzip_interleaved_process' in locals():
            if gzip_interleaved_process.wait() != 0:
                logging.error("gzip process for interleaved file failed.")
                sys.exit(1)
            logging.info("    Interleaved finished.")
            
    # Check if the output files were produced as expected
    check_file_exists(args.out_fastq_r1, "Output FASTQ file")
    check_file_exists(args.out_fastq_r2, "Output FASTQ file")
    check_file_exists(args.file_interleaved_gz, "Interleaved FASTQ file")

    if args.validate:
        logging.info("Validating FASTQ files...")
        validate_fastq(fastq_files_orig = [args.in_fastq_r1, args.in_fastq_r2], fastq_files_re = [args.out_fastq_r1, args.out_fastq_r2])

    if args.plot:
        # NB: order is important.  Plot script will use file names to help, otherwise it will assume order is R1, R2, R1-reoriented, R2-reoriented.
        fastq_files = [args.in_fastq_r1, args.in_fastq_r2, args.out_fastq_r1, args.out_fastq_r2]

        if args.plot_search_sequence:
            sequences = args.plot_search_sequence.upper().split(',')
        elif search_sequences_r1:
            # use all of the r1 sequences
            sequences = search_sequences_r1
        else:
            logging.error("Error: no search sequence available.  Either set --plot-search-sequence or --sequences-r1.")
            sys.exit(1)
        if args.plot_prefix:
            plot_prefix = args.plot_prefix
        else:
            plot_prefix = f"{output_dir}/result"
        plt.plot_positions(fastq_files, sequences = sequences, output_file_prefix=plot_prefix, n=100000)
    logging.info("Finished.")    

if __name__ == "__main__":
    main()
