import statistics
from collections import Counter
import argparse

# process phred scores in a fastq file and print statistics and top frequencies of scores

parser = argparse.ArgumentParser(description="Process a Phred string.")
parser.add_argument("phred_string", type=str, help="The Phred string to be processed")
args = parser.parse_args()

phred_string = args.phred_string

def phred_to_scores(phred_string):
    return [ord(char) - 33 for char in phred_string]

def calculate_statistics(scores):
    total = sum(scores)
    mean = statistics.mean(scores)
    median = statistics.median(scores)
    mode = statistics.mode(scores)
    
    return {
        'total': total,
        'mean': mean,
        'median': median,
        'mode': mode
    }

def get_frequencies(scores):
    counter = Counter(scores)
    most_common = counter.most_common(10)
    return most_common

scores = phred_to_scores(phred_string)
stats = calculate_statistics(scores)
frequencies = get_frequencies(scores)

print("Statistics:")
for key, value in stats.items():
    print(f"{key}: {value}")

print("\nScore frequencies (top 10):")
for value, count in frequencies:
    print(f"Score: {value}, Count: {count}")
