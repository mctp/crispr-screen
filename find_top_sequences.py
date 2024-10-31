import random
import subprocess
from Bio import SeqIO
from Bio import Align
from Bio.Align import AlignInfo
from Bio import motifs
import argparse
import gzip
import shutil
from difflib import SequenceMatcher
from termcolor import colored

parser = argparse.ArgumentParser(description="Identify top sequences in fastq and generate a consensus.")
parser.add_argument("--in-fastq", type=str, help="Input FASTQ file")
parser.add_argument("--out-fasta", type=str, help="Output FASTA file")
parser.add_argument("--out-alignment-clustalw", type=str, help="Output alignment file for ClustalW")
parser.add_argument("--sample-size", type=int, default=10, help="Number of sequences to sample")
parser.add_argument("--skip-consensus", action="store_true", help="Skip ClustalW alignment and consensus generation. Use if ClustalW is not available.")
args = parser.parse_args()

def is_gzipped(file_path):
    with open(file_path, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def sample_sequences(fastq_file, sample_size=10):
    sampled_fastq = "sampled.fastq"
    command = ["seqtk", "sample", fastq_file, str(sample_size)]
    
    if is_gzipped(fastq_file):
        command.insert(2, "-z")  # Indicate that the input is gzipped

    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"seqtk failed: {result.stderr}")

    with open(sampled_fastq, "w") as f:
        f.write(result.stdout)

    with open(sampled_fastq, "rt") as handle:
        sampled_sequences = list(SeqIO.parse(handle, "fastq"))
    
    return sampled_sequences

def write_fasta(sequences, output_file):
    with open(output_file, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")

def clustalw_alignment(fasta_file, alignment_file):
    command = ["clustalw2", "-infile=" + fasta_file, "-outfile=" + alignment_file, "-gapopen=100", "-gapext=10"]
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"ClustalW2 failed: {result.stderr}")
    alignment = Align.read(alignment_file, "clustal")
    return alignment

def get_consensus(alignment):
    motif = motifs.create(alignment)
    consensus = motif.degenerate_consensus
    return consensus

def determine_sequence_classes(sequences, threshold=0.9, min_fraction=0.1):
    sequence_classes = []
    total_sequences = len(sequences)
    
    for seq in sequences:
        added = False
        for seq_class in sequence_classes:
            if SequenceMatcher(None, seq_class[0].seq, seq.seq).ratio() >= threshold:
                seq_class.append(seq)
                added = True
                break
        if not added:
            sequence_classes.append([seq])
    
    # Prune classes with low representation
    pruned_classes = [seq_class for seq_class in sequence_classes if len(seq_class) / total_sequences >= min_fraction]
    
    return pruned_classes

sampled_sequences = sample_sequences(args.in_fastq, sample_size=args.sample_size)
write_fasta(sampled_sequences, args.out_fasta)

sequence_classes = determine_sequence_classes(sampled_sequences)

print(f"Number of sequence classes: {len(sequence_classes)}")
for i, seq_class in enumerate(sequence_classes):
    print(f"Class {i+1}:")
    print(f"Representative sequence: {seq_class[0].seq}")
    print(f"Number of sequences in class: {len(seq_class)} ({len(seq_class) / len(sampled_sequences):.2%})")

if not args.skip_consensus:
    clustalw_alignment = clustalw_alignment(args.out_fasta, args.out_alignment_clustalw)
    consensus = get_consensus(clustalw_alignment)

    if args.sample_size <= 10:
        print("Sampled sequences:")
        for sequence in sampled_sequences:
            print(sequence.seq)

        print("ClustalW2 Alignment:")
        print(clustalw_alignment)

    print("Consensus sequence:")
    def color_consensus(consensus):
        colored_consensus = ""
        for base in consensus:
            if base == 'A' or base == 'C' or base == 'T' or base == 'G':
                colored_consensus += colored(base, 'yellow', attrs=['underline'])
            elif base == 'N':
                colored_consensus += colored(base, 'red')
            else:
                colored_consensus += colored(base, 'blue')
        return colored_consensus

    print(color_consensus(consensus))
else:
    print("Consensus generation skipped.")
