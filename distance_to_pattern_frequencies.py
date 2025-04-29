import gzip
import sys
import argparse
from Bio import SeqIO
import subprocess
import logging
import shutil
# project-specific imports
import utilities.utilities as util

def is_available(name):
    try:
        result = subprocess.run([name, '--version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.returncode == 0
    except FileNotFoundError:
        return False

def open_fastq(file_path):
    if util.is_gzipped(file_path):
        if is_available('pigz'):
            return subprocess.Popen(['pigz', '-dc', file_path], stdout=subprocess.PIPE, text=True).stdout
        else:
            return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'rt')

def distance_to_pattern_frequencies(in_fastq=None, pattern=None, nreads=None, start_at=0, stop_at=10, frac_threshold=.01, show_full_sequence=False, center_sequence=False):  

    if frac_threshold is not None:
        logging.info(f"Using pruning threshold of {frac_threshold}")
    else:
        logging.info("Using default pruning threshold of 0.01")
        frac_threshold = .01

    if not all(c in 'ACGT' for c in pattern):
        logging.error("Pattern contains invalid characters. Only A, C, G, and T are allowed.")
        sys.exit(1)
        
    counter = 0
    distance_counts = {}
    representative_sequences = {}
    total_reads = 0
    with open_fastq(in_fastq) as handle:
        for i, record in enumerate(SeqIO.parse(handle, "fastq")):
            if nreads is not None and nreads != 0 and i >= nreads:
                break
            total_reads += 1
            sequence = str(record.seq)
            found = False
            tries = len(sequence) - len(pattern) + 1 if stop_at == 0 else min(stop_at, len(sequence) - len(pattern) + 1)
            for j in range(start_at, tries):
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
        logging.info(f"Pattern found in {counter} ({counter / total_reads * 100:.2f}%) of total reads.")
    sorted_distance_counts = dict(sorted(distance_counts.items(), key=lambda item: item[1], reverse=True))

    # lengths are for formatting purposes
    max_length = stop_at + len(pattern) + 25
    max_distance_length = max(len(f"Distance {distance}") for distance in sorted_distance_counts.keys())
    max_count_length = max(len(f"{count}") for count in sorted_distance_counts.values())
    max_percentage_length = max(len(f"{(count / total_reads) * 100:.2f}%") for count in sorted_distance_counts.values())

    terminal_width = shutil.get_terminal_size().columns
    logging.info(f"Terminal width: {terminal_width}")

    left_string_length = max_distance_length + 2 + max_count_length + 1 + max_percentage_length + 5
    center_pos = (terminal_width - left_string_length) // 2
    if center_sequence:
        max_distance = max(sorted_distance_counts.keys())
        if max_distance < center_pos:
            center_pos = max_distance

    logging.info("Distance frequencies and representative sequences:")
    for distance, count in sorted_distance_counts.items():
        fraction = count / total_reads
        percentage = fraction * 100
        if fraction <= frac_threshold:
            continue
        sequence = representative_sequences[distance]
        if distance == -1:
            distance_str = "Distance *"
        else:
            distance_str = f"Distance {distance}"
        
        distance_str = distance_str.ljust(max_distance_length)
        count_str = f"{count}".ljust(max_count_length)
        percentage_str = f"{percentage:.2f}%".ljust(max_percentage_length)
        log_prefix_length = len("[2025-03-24 20:00:24] [INFO] ") # example logging prefix length
        adjusted_terminal_width = terminal_width - log_prefix_length

        if distance == -1:
            if show_full_sequence:
                logging.info(f"{distance_str}: {count_str} ({percentage_str}) {sequence:<{max_length + 3}}")
            else:
                if len(sequence) > adjusted_terminal_width - left_string_length:
                    sequence = sequence[:adjusted_terminal_width - left_string_length - 3] + "..."
                logging.info(f"{distance_str}: {count_str} ({percentage_str}) {sequence}")
        else:
            if center_sequence:
                if distance < center_pos:
                    padding_left = center_pos - distance
                    sequence = " " * padding_left + sequence
                if center_pos + len(pattern) > len(sequence):
                    pattern_end = len(sequence) - center_pos
                    sequence = sequence[:center_pos] + f"\033[91m{pattern[:pattern_end]}\033[0m" + sequence[center_pos + pattern_end:]
                else:
                    sequence = sequence[:center_pos] + f"\033[91m{pattern}\033[0m" + sequence[center_pos + len(pattern):]
            else:
                sequence = sequence[:distance] + f"\033[91m{pattern}\033[0m" + sequence[distance + len(pattern):]
            if not show_full_sequence:
                if len(sequence) > adjusted_terminal_width - left_string_length:
                    if sequence.endswith(f"\033[91m{pattern}\033[0m"):
                        sequence = sequence[:adjusted_terminal_width - left_string_length + 6] + "\033[0m..."
                    else:
                        sequence = sequence[:adjusted_terminal_width - left_string_length + 6] + "..."
            logging.info(f"{distance_str}: {count_str} ({percentage_str}) {sequence}")

if __name__ == "__main__":
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
    parser = argparse.ArgumentParser(description="Filter sequences from a FASTQ file.")
    parser.add_argument('--in-fastq', required=True, help="Input FASTQ file (gzipped or plain text).")
    parser.add_argument('--pattern', required=True, default="GATCTAGTTACGCC", help="Pattern to search for in sequences.")
    parser.add_argument('--nreads', type=int, default=None, help="Number of reads in the FASTQ file to scan (default: all).")
    parser.add_argument('--stop-at', '--max-tries', type=int, default=10, help="Pattern should be fully captured within this many bases into the sequences (default: 10).")
    parser.add_argument('--start-at','--min-tries', type=int, default=0, help="Minimum number of bases into the sequences to start searching for the pattern (default: 0).")
    parser.add_argument('--show-full-sequence', action='store_true', help="Show full sequence in output.")
    parser.add_argument('--show-threshold', type=float, default=None, help="Show only distances with frequencies above this threshold.")
    parser.add_argument('--center-sequence','--center-sequences','--center', action='store_true', help="Center the pattern in the sequence output.")
    args = parser.parse_args()

    distance_to_pattern_frequencies(
        in_fastq=args.in_fastq, 
        pattern=args.pattern, 
        nreads=args.nreads, 
        start_at=args.start_at, 
        stop_at=args.stop_at, 
        frac_threshold=args.show_threshold, 
        show_full_sequence=args.show_full_sequence, 
        center_sequence=args.center_sequence)
