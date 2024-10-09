import argparse
import pysam
import sys
import multiprocessing
import os

def parse_args():
    parser = argparse.ArgumentParser(description='Analyze BAM file.')
    parser.add_argument('--bam', required=True, help='Path to the input BAM file')
    parser.add_argument('--cpus', type=int, default=1, help='Number of CPUs to use (default: 1)')
    return parser.parse_args()

args = parse_args()

def process_sgrna(sgrna_id, bam_file):
    bam = pysam.AlignmentFile(bam_file, "rb")
    sequences = {}
    mapq_frequencies = {}
    read_count = 0

    for alignment in bam.fetch(sgrna_id):
        alignment_id = alignment.query_name
        map_quality = alignment.mapping_quality
        sequence = alignment.query_sequence

        if map_quality in mapq_frequencies:
            mapq_frequencies[map_quality] += 1
        else:
            mapq_frequencies[map_quality] = 1

        if sgrna_id not in sequences:
            sequences[sgrna_id] = {"count": 0, "alignments": {}, "mapq_frequencies": {}}

        sequences[sgrna_id]["count"] += 1
        sequences[sgrna_id]["alignments"][alignment_id] = {
            "map_quality": map_quality,
            "sequence": sequence
        }
        read_count += 1

    if sgrna_id not in sequences:
        sequences[sgrna_id] = {"count": 0, "alignments": {}, "mapq_frequencies": {}}
    if not mapq_frequencies:
        mapq_frequencies = {"-1": 0}

    sequences[sgrna_id]["mapq_frequencies"] = mapq_frequencies

    bam.close()
    return sequences, read_count

def analyze_bam(bam_file, cpus):
    bam = pysam.AlignmentFile(bam_file, "rb")
    sgrna_ids = bam.references
    bam.close()

    pool = multiprocessing.Pool(processes=cpus)
    results = pool.starmap(process_sgrna, [(sgrna_id, bam_file) for sgrna_id in sgrna_ids])
    pool.close()
    pool.join()

    total_read_count = 0
    base_name = os.path.splitext(os.path.basename(bam_file))[0]
    sgrna_bam_results_file = f"{base_name}_sgrna_bam_results.txt"
    sgrna_bam_summary_file = f"{base_name}_sgrna_bam_summary.txt"
    sgrna_bam_freq_file = f"{base_name}_sgrna_bam_freq.txt"

    first_write = True # Flag to determine if we are writing the header
    for sequences, read_count in results:
        total_read_count += read_count

        # Determine the write mode
        mode = "w" if first_write else "a"

        # Print the result to sgrna_bam_results.txt
        with open(sgrna_bam_results_file, mode) as output_file:
            if first_write:
                output_file.write("sgrna_id\talignment_id\talignment_mapq\talignment_sequence\tcount\n")
            for sgrna_id, data in sequences.items():
                for alignment_id, alignment_data in data["alignments"].items():
                    output_file.write(f"{sgrna_id}\t{alignment_id}\t{alignment_data['map_quality']}\t{alignment_data['sequence']}\t{data['count']}\n")

        # Print the result to sgrna_bam_summary.txt
        with open(sgrna_bam_summary_file, mode) as highest_mapq_file:
            if first_write:
                highest_mapq_file.write("sgrna_id\ttotal_alignments\tmapq_1\tmapq_2\tmapq_1_frac\tmapq_2_frac\n")
            for sgrna_id, data in sequences.items():
                mapq_frequencies = data["mapq_frequencies"]
                if len(mapq_frequencies) > 0:
                    sorted_mapq = sorted(mapq_frequencies.items(), key=lambda item: item[1], reverse=True)
                    highest_mapq, highest_freq = sorted_mapq[0]
                    second_highest_mapq, second_highest_freq = sorted_mapq[1] if len(sorted_mapq) > 1 else (None, 0)
                    total_alignments = sum(mapq_frequencies.values())
                    if total_alignments > 0:
                        highest_fraction = highest_freq / total_alignments
                        second_highest_fraction = second_highest_freq / total_alignments if second_highest_mapq is not None else 0
                    else:
                        highest_fraction = 0
                        second_highest_fraction = 0
                    highest_mapq_file.write(f"{sgrna_id}\t{total_alignments}\t{highest_mapq}\t{second_highest_mapq}\t{highest_fraction:.4f}\t{second_highest_fraction:.4f}\n")

        # Print the frequencies to sgrna_bam_freq.txt
        with open(sgrna_bam_freq_file, mode) as freq_file:
            if first_write:
                freq_file.write("sgrna_id\tmapq\tfrequency\n")
            for sgrna_id, data in sequences.items():
                for mapq, frequency in data["mapq_frequencies"].items():
                    freq_file.write(f"{sgrna_id}\t{mapq}\t{frequency}\n")

        first_write = False

    print(f"Total reads: {total_read_count}\n")
    print(f"Results: {sgrna_bam_results_file}")
    print(f"Summary: {sgrna_bam_summary_file}")
    print(f"Frequencies: {sgrna_bam_freq_file}")

analyze_bam(args.bam, args.cpus)
