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
import utilities.parse_config as pc
import shutil
import time
# project-specific imports
import utilities.utilities as util


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

def gzip_file(input_file):
    if not os.path.exists(input_file):
        logging.error(f"Input file {input_file} does not exist.")
        raise FileNotFoundError(f"Input file {input_file} does not exist.")
    temp_file = f"{input_file}.tmp"
    try:
        with open(input_file, 'rb') as f_in, open(temp_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        if not os.path.exists(temp_file):
            logging.error(f"Temporary file {temp_file} was not created.")
            raise FileNotFoundError(f"Temporary file {temp_file} was not created.")
        if pigz_available:
            gzip = "pigz"
        else:
            gzip = "gzip"
        gzip_command = f"{gzip} {temp_file}"
        subprocess.run(gzip_command, shell=True, check=True)
        shutil.move(f"{temp_file}.gz", f"{input_file}.gz")
    except Exception as e:
        logging.error(f"Failed to gzip file {input_file}: {e}")
        if os.path.exists(temp_file):
            os.remove(temp_file)
        raise

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
    
def main(sequences_r1, sequences_r2, cpus, in_fastq_r1, in_fastq_r2, out_fastq_r1, out_fastq_r2, plot, plot_prefix, plot_search_sequence, plot_only, validate, batch_size=1000, chunk_size=1000000, nreads=None, file_interleaved_gz=None, config=None):
    global pigz_available
    pigz_available = util.pigz_available()
    global is_docker
    is_docker = util.is_docker_env()

    if not in_fastq_r1 or not in_fastq_r2 or not out_fastq_r1 or not out_fastq_r2:
        logging.error("Error: Missing required FASTQ file arguments.")
        sys.exit(1)

    if config:
        config = pc.parse_config(config)
        if is_docker:
            output_dir = config['DOCKER_OUTPUT_DIR']
        else:
            output_dir = config['OUTPUT_DIR']
    else:
        output_dir = "output"

    search_sequences_r1 = sequences_r1.upper().split(',')
    search_sequences_r2 = sequences_r2.upper().split(',') if sequences_r2 else []

    if not plot and plot_only:
        logging.error("Error: --plot-only specified without --plot.")
        sys.exit(1)

    if not plot_only:
        ncpu = int(cpus)
        logging.info(f"Checking input FASTQ file {in_fastq_r1}")
        check_fastq_file(in_fastq_r1)
        logging.info(f"Checking input FASTQ file {in_fastq_r2}")
        check_fastq_file(in_fastq_r2)

        ncpu = int(cpus)

        if not search_sequences_r1 or search_sequences_r1 == ['']:
            raise ValueError("R1 sequences must not be empty.")
        if not search_sequences_r2:
            logging.warning("R2 sequences are empty (allowed).")

        logging.info(f"R1 search sequences: {search_sequences_r1}")
        logging.info(f"R2 search sequences: {search_sequences_r2}")

        if file_interleaved_gz:
            fastq_interleaved = re.sub(r'\.gz$', '', file_interleaved_gz)
        else:
            if "R1" in in_fastq_r1:
                fastq_interleaved = re.sub(r'R1', 'interleaved', in_fastq_r1)
                fastq_interleaved = re.sub(r'\.(fastq|fq)\.gz$', '.fq', fastq_interleaved)
            else:
                r1_prefix = re.sub(r'\.(fastq|fq)(\.gz)?$', '', in_fastq_r1)
                fastq_interleaved = f"{r1_prefix}_interleaved.fq"
            file_interleaved_gz = f"{fastq_interleaved}.gz"

        if file_interleaved_gz != f"{fastq_interleaved}.gz":
            logging.error("Error: file_interleaved_gz must be identical to fastq_interleaved with '.gz' appended.")
            sys.exit(1)

        out_r1_re_fastq = out_fastq_r1[:-3]
        out_r2_re_fastq = out_fastq_r2[:-3]

        if os.path.exists(out_fastq_r1) and os.path.getsize(out_fastq_r1) > 0 and \
           os.path.exists(out_fastq_r2) and os.path.getsize(out_fastq_r2) > 0:
            logging.info(f"Final output files {out_fastq_r1} and {out_fastq_r2} already exist and are not empty. Skipping processing.")
        else:
            if not (os.path.exists(fastq_interleaved) and os.path.getsize(fastq_interleaved) > 0):
                logging.info(f"Interleaving FASTQ files {in_fastq_r1} and {in_fastq_r2} (using seqtk) ...")
                interleave_seqtk(in_fastq_r1, in_fastq_r2, fastq_interleaved, nreads=nreads)
                logging.info(f"Interleaved FASTQ written to {fastq_interleaved}")
                gzip_file(fastq_interleaved)
                logging.info(f"Interleaved FASTQ file {fastq_interleaved} or {file_interleaved_gz} already exists and is not empty.")

            if not os.path.exists(file_interleaved_gz) or os.path.getsize(file_interleaved_gz) == 0:
                logging.info(f"Compressing interleaved FASTQ file to {file_interleaved_gz}")
                gzip_file(fastq_interleaved)

            chunk_count = 0
            while os.path.exists(f"{fastq_interleaved}.{chunk_count}"):
                chunk_count += 1

            if chunk_count == 0:
                chunk_count = split_interleaved_file_seqtk(fastq_interleaved, chunk_size=chunk_size)
            logging.info(f"{chunk_count} chunked files generated from {fastq_interleaved}")

            if not (os.path.exists(out_r1_re_fastq) and os.path.getsize(out_r1_re_fastq) > 0) or \
               not (os.path.exists(out_r2_re_fastq) and os.path.getsize(out_r2_re_fastq) > 0):
                chunk_count = 0
                while os.path.exists(f"{fastq_interleaved}.{chunk_count}"):
                    chunk_file = f"{fastq_interleaved}.{chunk_count}"
                    out_r1_chunk = f"{out_r1_re_fastq}.{chunk_count}"
                    out_r2_chunk = f"{out_r2_re_fastq}.{chunk_count}"

                    if not (os.path.exists(out_r1_chunk) and os.path.getsize(out_r1_chunk) > 0) or \
                       not (os.path.exists(out_r2_chunk) and os.path.getsize(out_r2_chunk) > 0):
                        total_stats = reorient_fastq(chunk_file, out_r1_chunk, out_r2_chunk, search_sequences_r1=search_sequences_r1, search_sequences_r2=search_sequences_r2, cpus=ncpu, batch_size=batch_size)
                        for key, value in total_stats.items():
                            logging.info(f"{key}: {value}")
                    chunk_count += 1

                with open(out_r1_re_fastq, 'w') as f1_out, open(out_r2_re_fastq, 'w') as f2_out:
                    for i in range(chunk_count):
                        with open(f"{out_r1_re_fastq}.{i}", 'r') as f1_in, open(f"{out_r2_re_fastq}.{i}", 'r') as f2_in:
                            f1_out.writelines(f1_in.readlines())
                            f2_out.writelines(f2_in.readlines())

                logging.info(f"Reoriented FASTQ files written to {out_r1_re_fastq} and {out_r2_re_fastq}")
                gzip_file(out_r1_re_fastq)
                logging.info(f"Reoriented FASTQ files {out_r1_re_fastq} and {out_r2_re_fastq} already exist and are not empty.")

            if not os.path.exists(out_fastq_r1) or os.path.getsize(out_fastq_r1) == 0:
                if os.path.exists(out_r2_re_fastq) and os.path.getsize(out_r2_re_fastq) > 0:
                    gzip_file(out_r2_re_fastq)
                else:
                    logging.error(f"Error: File {out_r2_re_fastq} does not exist or is empty.")
                    raise FileNotFoundError(f"File {out_r2_re_fastq} does not exist or is empty.")
                gzip_file(out_r1_re_fastq)

            if not os.path.exists(out_fastq_r2) or os.path.getsize(out_fastq_r2) == 0:
                logging.info(f"Compressing {out_r2_re_fastq} to {out_fastq_r2}")
                gzip_file(out_r2_re_fastq)

            for i in range(chunk_count):
                os.remove(f"{fastq_interleaved}.{i}")
                os.remove(f"{out_r1_re_fastq}.{i}")
                os.remove(f"{out_r2_re_fastq}.{i}")

            os.remove(fastq_interleaved)
            os.remove(out_r1_re_fastq)
            os.remove(out_r2_re_fastq)
    else:
        logging.info("Skipping read reorientation.")
        input_fastq_files = [in_fastq_r1, in_fastq_r2, out_fastq_r1, out_fastq_r2]
        for fastq_file in input_fastq_files:
            if not os.path.exists(fastq_file):
                logging.error(f"Error: Input FASTQ file {fastq_file} does not exist.")
                sys.exit(1)
            if os.path.getsize(fastq_file) == 0:
                logging.error(f"Error: Input FASTQ file {fastq_file} is empty.")
                sys.exit(1)

    if validate:
        logging.info("Validating FASTQ files...")
        validate_fastq(fastq_files_orig=[in_fastq_r1, in_fastq_r2], fastq_files_re=[out_fastq_r1, out_fastq_r2])

    if plot:
        fastq_files = [in_fastq_r1, in_fastq_r2, out_fastq_r1, out_fastq_r2]
        sequences = plot_search_sequence.upper().split(',') if plot_search_sequence else search_sequences_r1
        plot_prefix = plot_prefix or f"{output_dir}/result"
        plt.plot_positions(fastq_files, sequences=sequences, output_file_prefix=plot_prefix, n=100000)

    logging.info("Finished.")

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Reorient FASTQ files.')
    parser.add_argument('--cpus','-p', type=int, default=1, help='Number of CPUs to use (default: 1). Must be a positive integer.')

    parser.add_argument('--sequences-r1', '--seq-r1', type=str, required=True, help='Search sequences for R1 (comma-separated)')
    parser.add_argument('--sequences-r2', '--seq-r2', type=str, required=True, help='Search sequences for R2 (comma-separated)')

    parser.add_argument('--batch-size', '-b', type=int, default=1000, help='Batch size for processing (default: 1000).')
    parser.add_argument('--chunk-size', '-c', type=int, default=1000000, help='Chunk size for splitting interleaved file (default: 1000000).')
    parser.add_argument('--nreads', '-n', type=int, default=None, help='Number of reads to process (default: all).')

    parser.add_argument('--in-fastq-r1', '--in-r1', type=str, required=True, help='Input R1 FASTQ file (gzipped).')
    parser.add_argument('--in-fastq-r2', '--in-r2', type=str, required=True, help='Input R2 FASTQ file (gzipped).')
    parser.add_argument('--out-fastq-r1', '--out-r1', type=str, required=True, help='Output reoriented R1 FASTQ file (gzipped).')
    parser.add_argument('--out-fastq-r2', '--out-r2', type=str, required=True, help='Output reoriented R2 FASTQ file (gzipped).')

    parser.add_argument('--file-interleaved-gz', '--int', '--interleaved-fastq', type=str, default=None, help='Output interleaved FASTQ file (gzipped).')

    parser.add_argument('--plot', action='store_true', help='Plot sequence type proportions.')
    parser.add_argument('--plot-prefix', type=str, help='Prefix for plot files.')
    parser.add_argument('--plot-search-sequence', '--plot-search-sequences', '--plot-seq', type=str, help='Search sequence(s) for plotting.')
    parser.add_argument('--plot-only', action='store_true', help='Skip read reorientation and do plots.')

    parser.add_argument('--validate', action='store_true', help='Validate FASTQ files after processing.')

    args = parser.parse_args()

    main(sequences_r1=args.sequences_r1,
        sequences_r2=args.sequences_r2,
        cpus=args.cpus,
        nreads=args.nreads,
        batch_size=args.batch_size,
        chunk_size=args.chunk_size,
        in_fastq_r1=args.in_fastq_r1,
        in_fastq_r2=args.in_fastq_r2,
        out_fastq_r1=args.out_fastq_r1,
        out_fastq_r2=args.out_fastq_r2,
        file_interleaved_gz=args.file_interleaved_gz,
        plot=args.plot,
        plot_prefix=args.plot_prefix,
        plot_search_sequence=args.plot_search_sequence,
        plot_only=args.plot_only,
        validate=args.validate
        )
