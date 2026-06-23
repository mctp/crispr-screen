#!/usr/bin/env python3
"""
Reorient paired-end FASTQ reads based on search patterns for CRISPR screen processing.
Version 2 updates that improve performance:
- Data streaming via Python iterators (no intermediate files)
- Batch processing with multiprocessing pool
- Parse fastq with pysam

Author: Brian Magnuson
Email: bmagnuso@umich.edu
Date: 2026-02-16
"""

import argparse
import gzip
import logging
import os
import subprocess
import sys
import time
from multiprocessing import Pool, cpu_count
from typing import Dict, Iterator, List, Tuple

# pysam is required
try:
    import pysam
except ImportError:
    print("ERROR: pysam is required but not installed.", file=sys.stderr)
    print("\nAttempting to install pysam...", file=sys.stderr)
    try:
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pysam"])
        import pysam
        print("Successfully installed pysam.", file=sys.stderr)
    except Exception as e:
        print(f"\nFailed to install pysam: {e}", file=sys.stderr)
        print("\nPlease install pysam manually:", file=sys.stderr)
        sys.exit(1)

# configure logging
logging.basicConfig(
    level=logging.INFO,
    format="[{asctime}] [{levelname}] {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)


# define container class for PE reads
# ReadPair:
#   r1_id, r1_seq, r1_qual
#   r2_id, r2_seq, r2_qual
class ReadPair:
    """Container class for paired-end reads."""
    __slots__ = ['r1_id', 'r1_seq', 'r1_qual', 'r2_id', 'r2_seq', 'r2_qual']
    
    def __init__(self, r1_id, r1_seq, r1_qual, r2_id, r2_seq, r2_qual):
        self.r1_id = r1_id
        self.r1_seq = r1_seq
        self.r1_qual = r1_qual
        self.r2_id = r2_id
        self.r2_seq = r2_seq
        self.r2_qual = r2_qual
    
    def to_fastq_r1(self):
        """Format R1 as FASTQ string."""
        return f"@{self.r1_id}\n{self.r1_seq}\n+\n{self.r1_qual}\n"
    
    def to_fastq_r2(self):
        """Format R2 as FASTQ string."""
        return f"@{self.r2_id}\n{self.r2_seq}\n+\n{self.r2_qual}\n"
    
    def swap(self):
        """Swap R1 and R2."""
        self.r1_id, self.r2_id = self.r2_id, self.r1_id
        self.r1_seq, self.r2_seq = self.r2_seq, self.r1_seq
        self.r1_qual, self.r2_qual = self.r2_qual, self.r1_qual


def stream_paired_reads(r1_path: str, r2_path: str, max_reads: int = None) -> Iterator[ReadPair]:
    """
    Stream paired reads from gzipped FASTQ files using pysam.
    
    Args:
        r1_path: Path to R1 gzipped FASTQ file
        r2_path: Path to R2 gzipped FASTQ file
        max_reads: Maximum number of read pairs to process (None = all)
    
    Yields:
        ReadPair objects
    """
    r1_reader = pysam.FastxFile(r1_path)
    r2_reader = pysam.FastxFile(r2_path)
    
    try:
        for i, (read1, read2) in enumerate(zip(r1_reader, r2_reader)):
            if max_reads and i >= max_reads:
                break
            
            # Add /1 and /2 suffixes to name if not present
            r1_name = read1.name if read1.name.endswith('/1') else f"{read1.name}/1"
            r2_name = read2.name if read2.name.endswith('/2') else f"{read2.name}/2"
            
            # full header: name/1 [comment]
            r1_header = r1_name
            if read1.comment:
                r1_header = f"{r1_name} {read1.comment}"
            
            r2_header = r2_name
            if read2.comment:
                r2_header = f"{r2_name} {read2.comment}"
            
            yield ReadPair(
                r1_id=r1_header,
                r1_seq=read1.sequence,
                r1_qual=read1.quality,
                r2_id=r2_header,
                r2_seq=read2.sequence,
                r2_qual=read2.quality
            )
    finally:
        r1_reader.close()
        r2_reader.close()


def batch_iterator(iterator: Iterator, batch_size: int) -> Iterator[List]:
    """
    Group items from iterator into batches (reduce multiprocessing overhead).
   
    Args:
        iterator: Source iterator
        batch_size: Number of items per batch
    
    Yields:
        Lists of items
    """
    batch = []
    for item in iterator:
        batch.append(item)
        if len(batch) >= batch_size:
            yield batch
            batch = []
    if batch:
        yield batch


def reporting_thresholds() -> Iterator[int]:
    """
    Reports every 100K reads until 1M, then every 1M reads thereafter.
    This is a generator that maintains state across calls.
    
    Yields:
        Read counts at which to log progress (100000, 200000, ..., 1000000, 2000000, ...)
    """
    threshold = 100_000
    interval_small = 100_000
    interval_large = 1_000_000
    
    while True:
        yield threshold
        if threshold < 1_000_000:
            threshold += interval_small
        else:
            threshold += interval_large


def process_batch(args: Tuple[List[ReadPair], List[str], List[str]]) -> Tuple[List[ReadPair], Dict[str, int]]:
    """
    Process a batch of read pairs - reorient (swap R1/R2) based on pattern matching.
    Matching order of operations:
    1. If R1 has any provided R1 pattern, keep current orientation.
    2. Else if R2 has any provided R2 pattern, keep current orientation.
    3. Else if R2 has any provided R1 pattern, swap R1/R2.
    4. Else if R1 has any provided R2 pattern, swap R1/R2.
    5. Else, keep original orientation (no patterns found, no action taken).
    
    Args:
        args: Tuple of (batch, patterns_r1, patterns_r2)
    
    Returns:
        Tuple of (processed_pairs, statistics)
    """
    batch, patterns_r1, patterns_r2 = args
    
    results = []
    stats = {
        'r1_fwd': 0,
        'r1_rev': 0,
        'r2_fwd': 0,
        'r2_rev': 0,
        'not_found': 0
    }
    
    for pair in batch:
        # find pattern occurrences
        r1_has_r1_pattern = any(p in pair.r1_seq for p in patterns_r1)
        r2_has_r2_pattern = patterns_r2 and any(p in pair.r2_seq for p in patterns_r2)
        r2_has_r1_pattern = any(p in pair.r2_seq for p in patterns_r1)
        r1_has_r2_pattern = patterns_r2 and any(p in pair.r1_seq for p in patterns_r2)
        
        if r1_has_r1_pattern:
            results.append(pair)
            stats['r1_fwd'] += 1
        elif r2_has_r2_pattern:
            results.append(pair)
            stats['r2_rev'] += 1
        elif r2_has_r1_pattern:
            pair.swap()
            results.append(pair)
            stats['r1_rev'] += 1
        elif r1_has_r2_pattern:
            pair.swap()
            results.append(pair)
            stats['r2_fwd'] += 1
        else:
            # no pattern found, so keep original orientation
            results.append(pair)
            stats['not_found'] += 1
    
    return results, stats


def write_compressed_fastq(output_r1: str, output_r2: str, 
                          results_iterator: Iterator[ReadPair],
                          use_pigz: bool = True,
                          buffer_size: int = 10000):
    """
    Write read pairs directly to compressed output via pipe.

    Args:
        output_r1: Path for R1 output (include .gz extension)
        output_r2: Path for R2 output (include .gz extension)
        results_iterator: Iterator of ReadPair objects
        use_pigz: Use pigz instead of gzip
        buffer_size: Number of reads to accumulate before writing
    """
    # Choose compressor
    compressor = 'pigz' if use_pigz else 'gzip'
    
    # Start compression processes
    proc_r1 = subprocess.Popen(
        [compressor, '-c'],
        stdin=subprocess.PIPE,
        stdout=open(output_r1, 'wb'),
        stderr=subprocess.PIPE
    )
    
    proc_r2 = subprocess.Popen(
        [compressor, '-c'],
        stdin=subprocess.PIPE,
        stdout=open(output_r2, 'wb'),
        stderr=subprocess.PIPE
    )
    
    try:
        # buffer writes
        r1_buffer = []
        r2_buffer = []
        
        for pair in results_iterator:
            r1_buffer.append(pair.to_fastq_r1())
            r2_buffer.append(pair.to_fastq_r2())
            
            if len(r1_buffer) >= buffer_size:
                proc_r1.stdin.write(''.join(r1_buffer).encode())
                proc_r2.stdin.write(''.join(r2_buffer).encode())
                r1_buffer = []
                r2_buffer = []
        
        # write remaining
        if r1_buffer:
            proc_r1.stdin.write(''.join(r1_buffer).encode())
            proc_r2.stdin.write(''.join(r2_buffer).encode())
        
        # close and wait
        proc_r1.stdin.close()
        proc_r2.stdin.close()
        proc_r1.wait()
        proc_r2.wait()
        
        if proc_r1.returncode != 0:
            raise RuntimeError(f"R1 compression failed: {proc_r1.stderr.read().decode()}")
        if proc_r2.returncode != 0:
            raise RuntimeError(f"R2 compression failed: {proc_r2.stderr.read().decode()}")
            
    except Exception as e:
        # Clean up on error
        proc_r1.kill()
        proc_r2.kill()
        raise e


def validate_inputs(r1_path: str, r2_path: str):
    """Validate input files exist and are readable."""
    for path in [r1_path, r2_path]:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Input file not found: {path}")
        if os.path.getsize(path) == 0:
            raise ValueError(f"Input file is empty: {path}")
        # Try to open as gzip
        try:
            with gzip.open(path, 'rt') as f:
                f.read(1)
        except Exception as e:
            raise ValueError(f"Cannot read gzipped file {path}: {e}")


def check_pigz_available() -> bool:
    """Check if pigz is available."""
    try:
        subprocess.run(['pigz', '--version'], 
                      capture_output=True, 
                      check=True,
                      timeout=5)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError, subprocess.TimeoutExpired):
        return False


def main(sequences_r1, sequences_r2='', in_fastq_r1=None, in_fastq_r2=None,
         out_fastq_r1=None, out_fastq_r2=None, plot=False, plot_prefix=None,
         cpus=1, batch_size=5000, nreads=None, force_gzip=False, log_file=None,
         plot_only=False, validate=False):

    if log_file:
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter(
            "[{asctime}] [{levelname}] {message}",
            style="{",
            datefmt="%Y-%m-%d %H:%M:%S"
        ))
        logging.getLogger().addHandler(file_handler)

    logging.info("="*60)
    logging.info("FASTQ Read Reorientation v2 (parallel)")
    logging.info("="*60)

    patterns_r1 = [s.strip().upper() for s in sequences_r1.split(',') if s.strip()]
    patterns_r2 = [s.strip().upper() for s in sequences_r2.split(',') if s.strip()] if sequences_r2 else []

    if not patterns_r1:
        raise ValueError("R1 sequences must not be empty")

    logging.info(f"R1 search patterns: {patterns_r1}")
    logging.info(f"R2 search patterns: {patterns_r2 if patterns_r2 else '(none)'}")

    logging.info("Validating input files...")
    validate_inputs(in_fastq_r1, in_fastq_r2)
    logging.info(f"  R1: {in_fastq_r1} ({os.path.getsize(in_fastq_r1) / 1e9:.2f} GB)")
    logging.info(f"  R2: {in_fastq_r2} ({os.path.getsize(in_fastq_r2) / 1e9:.2f} GB)")

    use_pigz = check_pigz_available() and not force_gzip
    logging.info(f"Compression: {'pigz' if use_pigz else 'gzip'}")
    logging.info("Parser: pysam")

    num_cpus = min(cpus, cpu_count())
    logging.info(f"Configuration: batch_size={batch_size}, cpus={num_cpus}")
    logging.info(f"Read limit: {nreads if nreads else 'none (process all)'}")

    logging.info("="*60)
    logging.info("Starting processing...")
    start_time = time.time()

    pairs = stream_paired_reads(in_fastq_r1, in_fastq_r2, max_reads=nreads)
    batches = batch_iterator(pairs, batch_size)

    total_stats = {}
    all_results = []
    processed_reads = 0

    reporter = reporting_thresholds()
    next_report_threshold = next(reporter)

    with Pool(processes=num_cpus) as pool:
        batch_args = ((batch, patterns_r1, patterns_r2) for batch in batches)

        for results, stats in pool.imap(process_batch, batch_args, chunksize=1):
            all_results.extend(results)
            processed_reads += len(results)

            for key, value in stats.items():
                total_stats[key] = total_stats.get(key, 0) + value

            if processed_reads >= next_report_threshold:
                elapsed = time.time() - start_time
                rate = processed_reads / elapsed if elapsed > 0 else 0
                logging.info(f"Processed {processed_reads:,} read pairs "
                           f"({rate:,.0f} reads/sec)")
                next_report_threshold = next(reporter)

    processing_time = time.time() - start_time
    logging.info(f"Processing complete: {processed_reads:,} read pairs in {processing_time:.1f}s")
    logging.info(f"Rate: {processed_reads/processing_time:,.0f} read pairs/sec")

    logging.info("Reorientation statistics:")
    for key, value in sorted(total_stats.items()):
        pct = 100 * value / processed_reads if processed_reads > 0 else 0
        logging.info(f"  {key}: {value:,} ({pct:.1f}%)")

    logging.info("Writing compressed outputs...")
    write_start = time.time()
    write_compressed_fastq(
        out_fastq_r1,
        out_fastq_r2,
        iter(all_results),
        use_pigz=use_pigz
    )
    write_time = time.time() - write_start

    logging.info(f"Writing complete in {write_time:.1f}s")
    logging.info(f"  R1: {out_fastq_r1} ({os.path.getsize(out_fastq_r1) / 1e9:.2f} GB)")
    logging.info(f"  R2: {out_fastq_r2} ({os.path.getsize(out_fastq_r2) / 1e9:.2f} GB)")

    if plot:
        import read_reorientation_barplot as barplot
        fastq_files = [in_fastq_r1, in_fastq_r2, out_fastq_r1, out_fastq_r2]
        barplot.plot_positions(fastq_files, sequences=patterns_r1,
                               output_file_prefix=plot_prefix, n=100000)

    total_time = time.time() - start_time
    logging.info("="*60)
    logging.info(f"Total time: {total_time:.1f}s ({total_time/60:.1f} min)")
    logging.info(f"Overall rate: {processed_reads/total_time:,.0f} read pairs/sec")
    logging.info("DONE")

    return total_stats


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Reorient paired-end FASTQ reads based on search patterns (streaming version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--in-fastq-r1', '--in-r1', type=str, required=True,
                       help='Input R1 FASTQ file (gzipped)')
    parser.add_argument('--in-fastq-r2', '--in-r2', type=str, required=True,
                       help='Input R2 FASTQ file (gzipped)')
    parser.add_argument('--out-fastq-r1', '--out-r1', type=str, required=True,
                       help='Output reoriented R1 FASTQ file (gzipped)')
    parser.add_argument('--out-fastq-r2', '--out-r2', type=str, required=True,
                       help='Output reoriented R2 FASTQ file (gzipped)')
    parser.add_argument('--sequences-r1', '--seq-r1', type=str, required=True,
                       help='Search sequences for R1 (comma-separated)')
    parser.add_argument('--sequences-r2', '--seq-r2', type=str, default='',
                       help='Search sequences for R2 (comma-separated, optional)')
    parser.add_argument('--cpus', '-p', type=int, default=1,
                       help='Number of CPUs to use')
    parser.add_argument('--batch-size', '-b', type=int, default=5000,
                       help='Batch size for multiprocessing')
    parser.add_argument('--nreads', '-n', type=int, default=None,
                       help='Limit number of read pairs to process (for testing)')
    parser.add_argument('--force-gzip', action='store_true',
                       help='Use gzip instead of pigz')
    parser.add_argument('--log-file', type=str,
                       help='Write output to a file')
    parser.add_argument('--plot', action='store_true',
                       help='Plot sequence type proportions')
    parser.add_argument('--plot-prefix', type=str,
                       help='Prefix for plot files')

    args = parser.parse_args()

    if os.path.exists(args.out_fastq_r1) and os.path.exists(args.out_fastq_r2):
        if os.path.getsize(args.out_fastq_r1) > 0 and os.path.getsize(args.out_fastq_r2) > 0:
            logging.warning("Output files already exist and are not empty:")
            logging.warning(f"  {args.out_fastq_r1}")
            logging.warning(f"  {args.out_fastq_r2}")
            response = input("Overwrite? [y/N]: ")
            if response.lower() != 'y':
                logging.info("Exiting without processing.")
                sys.exit(0)

    try:
        main(
            sequences_r1=args.sequences_r1,
            sequences_r2=args.sequences_r2,
            in_fastq_r1=args.in_fastq_r1,
            in_fastq_r2=args.in_fastq_r2,
            out_fastq_r1=args.out_fastq_r1,
            out_fastq_r2=args.out_fastq_r2,
            plot=args.plot,
            plot_prefix=args.plot_prefix,
            cpus=args.cpus,
            batch_size=args.batch_size,
            nreads=args.nreads,
            force_gzip=args.force_gzip,
            log_file=args.log_file,
        )
    except KeyboardInterrupt:
        logging.error("Interrupted by user")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error: {e}", exc_info=True)
        sys.exit(1)
