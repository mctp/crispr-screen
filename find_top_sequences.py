import subprocess
from Bio import SeqIO
from Bio import Align
from Bio import motifs
import argparse
from difflib import SequenceMatcher
import os
import joblib
import numpy as np
import logging
import sys
# project-specific imports
import utilities.utilities as util

def sample_sequences(fastq_file, sample_size=10, seed=None):
    sampled_fastq = "sampled.fastq"
    command = ["seqtk", "sample", fastq_file, str(sample_size)]
    
    if util.is_gzipped(fastq_file):
        command.insert(2, "-z")

    if seed is not None:
        command.extend(["-s", str(seed)])

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
            if SequenceMatcher(None, seq_class[0][0].seq, seq.seq).ratio() >= threshold:
                seq_class.append((seq, None))
                added = True
                break
        if not added:
            sequence_classes.append([(seq, None)])
    
    # Prune classes with low representation
    pruned_classes = [seq_class for seq_class in sequence_classes if len(seq_class) / total_sequences >= min_fraction]
    
    return pruned_classes

def extract_conserved_kmers(sequence, k=6, var_start=65, var_end=85):
    left_kmers = [sequence[i:i+k] for i in range(var_start - k)]
    right_kmers = [sequence[i:i+k] for i in range(var_end, len(sequence) - k + 1)]
    return left_kmers + right_kmers

def create_feature_vector(sequence, top_kmers, var_start=65, var_end=85):
    kmers = extract_conserved_kmers(str(sequence.seq), var_start=var_start, var_end=var_end)
    return [kmers.count(kmer) for kmer in top_kmers]

# Classification phase
def load_model_and_classify(sequences, model_filename='dna_classifier.joblib', min_fraction=0.1):
    # Load the model and top_kmers
    loaded_data = joblib.load(model_filename)
    clf = loaded_data['model']
    if 'top_kmers' not in loaded_data:
        raise KeyError("The loaded model data does not contain 'top_kmers'")
    top_kmers = loaded_data['top_kmers']
    confidence_threshold = loaded_data['confidence_threshold']
    logging.info(f"Confidence threshold: {confidence_threshold}")

    # Create feature vectors for new sequences
    X = [create_feature_vector(seq, top_kmers) for seq in sequences]

    y_pred_proba = clf.predict_proba(X)
    
    # Get max probability for each prediction
    max_proba = np.max(y_pred_proba, axis=1)
    
    # Classify
    predictions = clf.predict(X)
    
    # Detect OOD samples
    is_ood = max_proba < confidence_threshold
    
    # Group sequences by predicted class
    sequence_classes = []
    for seq, pred, ood in zip(sequences, predictions, is_ood):
        added = False
        for seq_class in sequence_classes:
            if seq_class[0][1] == pred:
                seq_class.append((seq, pred, ood))
                added = True
                break
        if not added:
            sequence_classes.append([(seq, pred, ood)])

    return sequence_classes

def assign_representative_sequences(sequence_classes):
    representative_sequences = []
    for seq_class in sequence_classes:
        sequence_counts = {}
        for item in seq_class:
            seq = item[0]
            seq_str = str(seq.seq)
            if seq_str in sequence_counts:
                sequence_counts[seq_str] += 1
            else:
                sequence_counts[seq_str] = 1
        most_frequent_sequence = max(sequence_counts, key=sequence_counts.get)
        representative_sequences.append(most_frequent_sequence)
    return representative_sequences

def find_top_sequences(in_fastq, out_fasta, out_alignment_clustalw=None, sample_size=1000, msa_sample_size=200, skip_consensus=False, seed=None, method_2=False, model=None):
    try:
        from termcolor import colored
        termcolor_available = True
    except ImportError:
        termcolor_available = False

    base_name, ext1 = os.path.splitext(in_fastq)
    if ext1 == ".gz":
        base_name, ext2 = os.path.splitext(base_name)

    if not out_fasta:
        out_fasta = base_name + ".fasta"
    else:
        out_fasta = out_fasta

    if out_alignment_clustalw is None:
        out_alignment_clustalw = base_name + ".aln"

    sampled_sequences = sample_sequences(in_fastq, sample_size=sample_size, seed=seed)
    write_fasta(sampled_sequences, out_fasta)

    if method_2:
        if model:
            sequence_classes = load_model_and_classify(sampled_sequences, model_filename=model)
        else:
            sequence_classes = load_model_and_classify(sampled_sequences)
        ood_sequences = [seq for seq_class in sequence_classes for seq, _, ood in seq_class if ood]
        logging.info(f"Number of OOD sequences: {len(ood_sequences)} ({len(ood_sequences) / len(sampled_sequences):.2%})")
        filtered_sequence_classes = []
        for seq_class in sequence_classes:
            filtered_class = [seq for seq in seq_class if not seq[2]]
            if filtered_class:
                filtered_sequence_classes.append(filtered_class)
        
        logging.info("OOD sequences (first 20):")
        for ood_seq in ood_sequences[:20]:
            logging.info(ood_seq.seq)
        original_sequence_classes = sequence_classes
        sequence_classes = filtered_sequence_classes
    else:
        sequence_classes = determine_sequence_classes(sampled_sequences)
        original_sequence_classes = sequence_classes

    representative_sequences = assign_representative_sequences(sequence_classes)

    logging.info(f"Number of sequence classes: {len(sequence_classes)}")
    for i, seq_class in enumerate(sequence_classes):
        logging.info(f"Class {i+1}:")
        representative_sequence = representative_sequences[i]
        logging.info(f"Representative sequence: {representative_sequence}")
        logging.info(f"Number of sequences in class: {len(seq_class)} ({len(seq_class) / len(sampled_sequences):.2%})")


    if not skip_consensus:
        sampled_sequences_for_msa = sampled_sequences[:msa_sample_size]
        write_fasta(sampled_sequences_for_msa, out_fasta)
        clustalw_alignment_string = clustalw_alignment(out_fasta, out_alignment_clustalw)
        consensus = get_consensus(clustalw_alignment_string)

        if sample_size <= 10:
            logging.info("Sampled sequences:")
            for sequence in sampled_sequences:
                logging.info(sequence.seq)

            logging.info("ClustalW2 Alignment:")
            logging.info(clustalw_alignment_string)

        logging.info("Consensus sequence:")
        if termcolor_available:
            colored_consensus = ""
            for base in consensus:
                if base == 'A' or base == 'C' or base == 'T' or base == 'G':
                    colored_consensus += colored(base, 'yellow', attrs=['underline'])
                elif base == 'N':
                    colored_consensus += colored(base, 'red')
                else:
                    colored_consensus += colored(base, 'blue')
            logging.info(colored_consensus)
        else:
            logging.info(consensus)
    else:
        logging.info("Consensus generation skipped.")

    logging.info(f"FASTA file path: {out_fasta}")
    logging.info(f"Alignment file path: {out_alignment_clustalw}")


if __name__ == "__main__":

    logging.basicConfig(
        level=logging.INFO,
        format="[{asctime}] [{levelname}] {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler("output/find_top_sequences.log")
        ]
    )

    parser = argparse.ArgumentParser(description="Identify top sequences in fastq and generate a consensus.")
    parser.add_argument("--in-fastq", type=str, help="Input FASTQ file")
    parser.add_argument("--out-fasta", type=str, help="Output FASTA file")
    parser.add_argument("--out-alignment-clustalw", type=str, help="Output alignment file for ClustalW")
    parser.add_argument("--sample-size", type=int, default=1000, help="Number of sequences to sample")
    parser.add_argument("--msa-sample-size", type=int, default=200, help="Number of sequences to sample for MSA.")
    parser.add_argument("--skip-consensus", action="store_true", help="Skip ClustalW alignment and consensus generation. Use if ClustalW is not available.")
    parser.add_argument("--seed", type=int, help="Seed for random sampling")
    parser.add_argument("--method-2", action="store_true", help="Use the second method to classify sequences")
    parser.add_argument('--model', help='Filename of trained model.')
    args = parser.parse_args()

    find_top_sequences(
        in_fastq=args.in_fastq,
        out_fasta=args.out_fasta,
        out_alignment_clustalw=args.out_alignment_clustalw,
        sample_size=args.sample_size,
        msa_sample_size=args.msa_sample_size,
        skip_consensus=args.skip_consensus,
        seed=args.seed,
        method_2=args.method_2,
        model=args.model
   )
