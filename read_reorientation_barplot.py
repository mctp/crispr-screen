from Bio import SeqIO
import matplotlib.pyplot as plt
import gzip
import logging
import sys

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
def find_sequence_positions(fastq_files, sequences, n=None):
    data = {}
    for file in fastq_files:
        with gzip.open(file, "rt") as handle:
            if "R1" in file and "reoriented" not in file:
                id = "R1 original"
            elif "R1" in file and "reoriented" in file:
                id = "R1 reoriented"
            elif "R2" in file and "reoriented" not in file:
                id = "R2 original"
            elif "R2" in file and "reoriented" in file:
                id = "R2 reoriented"
            else:
                id = "R1 original" if "R1" in file else "R2 original"

            if file not in data:
                data[file] = {}
            data[file]['id'] = id
            data[file]['data'] = {}
            data[file]['any_seq_hits'] = 0
            data[file]['total'] = 0

            logging.info(f"Searching for sequence(s) {sequences} in {file}.")

            for record in SeqIO.parse(handle, "fastq"):
                data[file]['total'] += 1
                for seq in sequences:
                    if seq not in data[file]['data']:
                        data[file]['data'][seq] = {}
                        data[file]['data'][seq]['hits'] = 0
                    pos = str(record.seq).find(seq)
                    if pos not in data[file]['data'][seq]:
                        data[file]['data'][seq][pos] = 0
                    data[file]['data'][seq][pos] += 1
                    if pos != -1:
                        data[file]['data'][seq]['hits'] += 1
                        data[file]['any_seq_hits'] += 1
                        break
                if n and data[file]['total'] >= n:
                    break
    return data

def plot_positions(fastq_files, sequences, output_file_prefix="plot", n=100000):
    results_file = f"{output_file_prefix}_results.tsv"
    plot_png_file = f"{output_file_prefix}_barplot.png"
    plot_pdf_file = f"{output_file_prefix}_barplot.pdf"

    logging.info(f"Analyzing {n} reads per fastq file.")

    positions = find_sequence_positions(fastq_files, sequences, n=n) 
    
    # save results to
    with open(results_file, 'w') as f:
        f.write("File\tID\tSequence\tPosition\tCount\tSeqHits\tTotal Sequences\tAnySeqHits\n")
        for file, data in positions.items():
            id = data['id']
            total = data['total']
            any_seq_hits = data['any_seq_hits']
            data = data['data']
            for seq, seqdata in data.items():
                seq_hits = seqdata['hits']
                seqdata = dict(sorted(seqdata.items(), key=lambda item: item[1], reverse=True))
                for pos, count in seqdata.items():
                    if pos == 'hits':
                        continue
                    f.write(f"{file}\t{id}\t{seq}\t{pos}\t{count}\t{seq_hits}\t{total}\t{any_seq_hits}\n")
    logging.info(f"Results saved to {results_file}")
    # init plot
    _, ax = plt.subplots(figsize=(10, 6))

    classes = ['Found', 'Not Found']
    colors = plt.cm.get_cmap('tab20', len(classes))
    color_map = {cls: colors(i) for i, cls in enumerate(classes)}

    for file, data in positions.items():
        id = data['id']
        total_sequences = data['total']
        found_sequences = data['any_seq_hits']

        not_found_sequences = total_sequences - found_sequences

        found_percentage = (found_sequences / total_sequences) * 100
        not_found_percentage = (not_found_sequences / total_sequences) * 100

        percentages = [found_percentage, not_found_percentage]
        labels = ['Found', 'Not Found']

        logging.info(f"{id}:")
        for label, perc in zip(labels, percentages):
            logging.info(f"{label}: {perc:.2f}%")

        primary_sequence = True
        for seq, seqdata in data['data'].items():
            seq_hits = seqdata['hits']
            seq_hits_percentage = (seq_hits / total_sequences) * 100
            if primary_sequence:
                logging.info(f"Primary Sequence: {seq}, Hits: {seq_hits}, Percentage: {seq_hits_percentage:.2f}%")
                primary_sequence = False
            else:
                logging.info(f"Secondary Sequence: {seq}, Hits: {seq_hits}, Percentage: {seq_hits_percentage:.2f}%")

        bottom = 0
        for label, perc in zip(labels, percentages):
            bar = ax.bar(id, perc, bottom=bottom, color=color_map[label], label=label)
            if label == 'Found':
                ax.text(bar[0].get_x() + bar[0].get_width() / 2, bottom + perc / 2, f'{perc:.2f}%', ha='center', va='center', color='white', fontsize=10, fontweight='bold')
            bottom += perc

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    unique_labels = {label: handle for label, handle in by_label.items()}
    ax.legend(unique_labels.values(), unique_labels.keys(), title='Class', bbox_to_anchor=(1.05, 1), loc='upper left')

    ax.set_xlabel('Fastq file')
    ax.set_ylabel('Percentage')
    ax.set_title(f'Reads containing any of {sequences}')
    plt.tight_layout()
    plt.savefig(plot_png_file)
    logging.info(f"Plot saved to {plot_png_file}")
    plt.savefig(plot_pdf_file)
    logging.info(f"Plot saved to {plot_pdf_file}")
