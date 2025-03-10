import argparse
import os
import joblib
from collections import Counter
import numpy as np
from sklearn.naive_bayes import MultinomialNB
from sklearn.model_selection import train_test_split

parser = argparse.ArgumentParser(description='Train a DNA sequence classifier.')
parser.add_argument('--input', required=True, help='Path to the input file containing DNA sequences and their classes.')
parser.add_argument('--model-name', required=True, help='Filename to save the trained model.')
args = parser.parse_args()

# Validate input file
if not os.path.isfile(args.input):
    raise FileNotFoundError(f"Input file {args.input} does not exist.")

# Validate input file format
with open(args.input, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) != 2:
            raise ValueError(f"Input file {args.input} is not in the expected format. Each line must have exactly 2 tab-delimited columns.")

def load_data(input_file):
    sequences = []
    labels = []
    with open(input_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            sequences.append(parts[0])
            labels.append(int(parts[1]))
    return sequences, labels

def extract_conserved_kmers(sequence, k=6, var_start=65, var_end=85):
    left_kmers = [sequence[i:i+k] for i in range(var_start - k)]
    right_kmers = [sequence[i:i+k] for i in range(var_end, len(sequence) - k + 1)]
    return left_kmers + right_kmers

def create_feature_vector(sequence, top_kmers, var_start=65, var_end=85):
    kmers = extract_conserved_kmers(sequence, var_start=var_start, var_end=var_end)
    return [kmers.count(kmer) for kmer in top_kmers]

# Training phase
def train_and_save_model(sequences, labels, model_filename='dna_classifier.joblib'):
    kmers = [extract_conserved_kmers(seq) for seq in sequences]
    all_kmers = [kmer for seq_kmers in kmers for kmer in seq_kmers]
    top_kmers = [k for k, _ in Counter(all_kmers).most_common(100)]

    X = [create_feature_vector(seq, top_kmers) for seq in sequences]
    X_train, X_test, y_train, y_test = train_test_split(X, labels, test_size=0.2)

    clf = MultinomialNB()
    clf.fit(X_train, y_train)

    accuracy = clf.score(X_test, y_test)
    print(f"Model accuracy: {accuracy}")

# Calculate confidence threshold
    y_pred_proba = clf.predict_proba(X_test)
    confidence_threshold = np.percentile(np.max(y_pred_proba, axis=1), 5)
    
    # Save the model, top_kmers, and confidence threshold
    joblib.dump({'model': clf, 'confidence_threshold': confidence_threshold, 'top_kmers': top_kmers}, model_filename)

    print(f"Model saved to {model_filename}")

sequences, labels = load_data(args.input)
train_and_save_model(sequences, labels, model_filename=args.model_name)

