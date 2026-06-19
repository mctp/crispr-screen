from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.cm
import gzip
import logging
import sys
import argparse
import os
import glob

# Import project-specific modules for config parsing (only when needed)
try:
    import utilities.parse_config as pc
    import utilities.process_metadata as pm
except ImportError:
    # Fallback if utilities not available - script can still work with explicit parameters
    pc = None
    pm = None

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

def find_sequence_positions_separate(fastq_files, r1_sequences, r2_sequences, n=None):
    """Find sequence positions with separate tracking of R1 and R2 sequences"""
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
            data[file]['r1_hits'] = 0
            data[file]['r2_hits'] = 0
            data[file]['total'] = 0

            logging.info(f"Searching for R1 sequences {r1_sequences} and R2 sequences {r2_sequences} in {file}.")

            for record in SeqIO.parse(handle, "fastq"):
                data[file]['total'] += 1
                found_in_read = False
                
                # Check R1 sequences first
                for seq in r1_sequences:
                    if seq not in data[file]['data']:
                        data[file]['data'][seq] = {}
                        data[file]['data'][seq]['hits'] = 0
                    pos = str(record.seq).find(seq)
                    if pos not in data[file]['data'][seq]:
                        data[file]['data'][seq][pos] = 0
                    data[file]['data'][seq][pos] += 1
                    if pos != -1:
                        data[file]['data'][seq]['hits'] += 1
                        if not found_in_read:
                            data[file]['r1_hits'] += 1
                            data[file]['any_seq_hits'] += 1
                            found_in_read = True
                        break  # Only count first match per read
                
                # If no R1 sequence found, check R2 sequences
                if not found_in_read:
                    for seq in r2_sequences:
                        if seq not in data[file]['data']:
                            data[file]['data'][seq] = {}
                            data[file]['data'][seq]['hits'] = 0
                        pos = str(record.seq).find(seq)
                        if pos not in data[file]['data'][seq]:
                            data[file]['data'][seq][pos] = 0
                        data[file]['data'][seq][pos] += 1
                        if pos != -1:
                            data[file]['data'][seq]['hits'] += 1
                            data[file]['r2_hits'] += 1
                            data[file]['any_seq_hits'] += 1
                            break  # Only count first match per read
                
                if n and data[file]['total'] >= n:
                    break
    return data

def plot_positions(fastq_files, sequences, output_file_prefix="plot", n=100000, sample_name=None, library_name=None, r1_sequences=None, r2_sequences=None):
    results_file = f"{output_file_prefix}_results.tsv"
    plot_png_file = f"{output_file_prefix}_barplot.png"
    plot_pdf_file = f"{output_file_prefix}_barplot.pdf"

    logging.info(f"Analyzing {n} reads per fastq file.")

    # Determine if we're using separate R1/R2 sequences or combined sequences
    use_separate_r1_r2 = r1_sequences is not None and r2_sequences is not None
    
    if use_separate_r1_r2:
        positions = find_sequence_positions_separate(fastq_files, r1_sequences, r2_sequences, n=n)
        plot_classes = ['R1 Found', 'R2 Found', 'Not Found']
        logging.info(f"Using separate R1 sequences: {r1_sequences}")
        logging.info(f"Using separate R2 sequences: {r2_sequences}")
    else:
        positions = find_sequence_positions(fastq_files, sequences, n=n)
        plot_classes = ['Found', 'Not Found']
        logging.info(f"Using combined sequences: {sequences}")
    
    # save results to
    with open(results_file, 'w') as f:
        if use_separate_r1_r2:
            f.write("Library\tSample\tFile\tID\tSequence_Type\tSequence\tPosition\tCount\tSeqHits\tTotal Sequences\tR1_Hits\tR2_Hits\tAnySeqHits\n")
        else:
            f.write("Library\tSample\tFile\tID\tSequence\tPosition\tCount\tSeqHits\tTotal Sequences\tAnySeqHits\n")
        
        for file, data in positions.items():
            id = data['id']
            total = data['total']
            lib_col = library_name if library_name else "Unknown"
            sample_col = sample_name if sample_name else "Unknown"
            
            if use_separate_r1_r2:
                # Write separate R1 and R2 data
                any_seq_hits = data['any_seq_hits']
                r1_hits = data.get('r1_hits', 0)
                r2_hits = data.get('r2_hits', 0)
                
                # Write R1 sequence data
                for seq in r1_sequences:
                    if seq in data['data']:
                        seqdata = data['data'][seq]
                        seq_hits = seqdata['hits']
                        seqdata_sorted = dict(sorted(seqdata.items(), key=lambda item: item[1], reverse=True))
                        for pos, count in seqdata_sorted.items():
                            if pos == 'hits':
                                continue
                            f.write(f"{lib_col}\t{sample_col}\t{file}\t{id}\tR1\t{seq}\t{pos}\t{count}\t{seq_hits}\t{total}\t{r1_hits}\t{r2_hits}\t{any_seq_hits}\n")
                
                # Write R2 sequence data
                for seq in r2_sequences:
                    if seq in data['data']:
                        seqdata = data['data'][seq]
                        seq_hits = seqdata['hits']
                        seqdata_sorted = dict(sorted(seqdata.items(), key=lambda item: item[1], reverse=True))
                        for pos, count in seqdata_sorted.items():
                            if pos == 'hits':
                                continue
                            f.write(f"{lib_col}\t{sample_col}\t{file}\t{id}\tR2\t{seq}\t{pos}\t{count}\t{seq_hits}\t{total}\t{r1_hits}\t{r2_hits}\t{any_seq_hits}\n")
            else:
                # Original format for backward compatibility
                any_seq_hits = data['any_seq_hits']
                data_dict = data['data']
                for seq, seqdata in data_dict.items():
                    seq_hits = seqdata['hits']
                    seqdata_sorted = dict(sorted(seqdata.items(), key=lambda item: item[1], reverse=True))
                    for pos, count in seqdata_sorted.items():
                        if pos == 'hits':
                            continue
                        f.write(f"{lib_col}\t{sample_col}\t{file}\t{id}\t{seq}\t{pos}\t{count}\t{seq_hits}\t{total}\t{any_seq_hits}\n")
    
    logging.info(f"Results saved to {results_file}")
    
    # init plot
    _, ax = plt.subplots(figsize=(12, 6))

    if use_separate_r1_r2:
        # Custom colors for R1/R2 mode: blue for R1, light blue for R2, light red for not found
        color_map = {
            'R1 Found': '#1f77b4',     # Blue
            'R2 Found': '#87ceeb',     # Light blue (sky blue)
            'Not Found': '#ffcccb'     # Light red
        }
    else:
        # Original color scheme for backward compatibility
        colors = matplotlib.cm.get_cmap('tab20')
        color_map = {cls: colors(i) for i, cls in enumerate(plot_classes)}

    # Use positional indices for proper bar spacing
    file_ids = []
    x_positions = []
    for idx, (file, data) in enumerate(positions.items()):
        id = data['id']
        file_ids.append(id)
        x_positions.append(idx)
        total_sequences = data['total']
        
        if use_separate_r1_r2:
            r1_hits = data.get('r1_hits', 0)
            r2_hits = data.get('r2_hits', 0)
            not_found = total_sequences - r1_hits - r2_hits
            
            r1_percentage = (r1_hits / total_sequences) * 100 if total_sequences > 0 else 0
            r2_percentage = (r2_hits / total_sequences) * 100 if total_sequences > 0 else 0
            not_found_percentage = (not_found / total_sequences) * 100 if total_sequences > 0 else 0
            
            percentages = [r1_percentage, r2_percentage, not_found_percentage]
            labels = ['R1 Found', 'R2 Found', 'Not Found']
            
            logging.info(f"{id}:")
            logging.info(f"R1 Found: {r1_percentage:.2f}%")
            logging.info(f"R2 Found: {r2_percentage:.2f}%")
            logging.info(f"Not Found: {not_found_percentage:.2f}%")
            
            # Log individual R1 sequences
            for seq in r1_sequences:
                if seq in data['data']:
                    seq_hits = data['data'][seq]['hits']
                    seq_percentage = (seq_hits / total_sequences) * 100 if total_sequences > 0 else 0
                    logging.info(f"R1 Sequence {seq}: {seq_hits} hits ({seq_percentage:.2f}%)")
            
            # Log individual R2 sequences
            for seq in r2_sequences:
                if seq in data['data']:
                    seq_hits = data['data'][seq]['hits']
                    seq_percentage = (seq_hits / total_sequences) * 100 if total_sequences > 0 else 0
                    logging.info(f"R2 Sequence {seq}: {seq_hits} hits ({seq_percentage:.2f}%)")
        else:
            # Original two-class approach
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
            bar = ax.bar(idx, perc, bottom=bottom, color=color_map[label], label=label, width=0.6)
            if perc > 5:  # Only show text if percentage is > 5%
                ax.text(bar[0].get_x() + bar[0].get_width() / 2, bottom + perc / 2, f'{perc:.1f}%', 
                       ha='center', va='center', color='white', fontsize=9, fontweight='bold')
            bottom += perc
    
    # Set x-axis to use positional indices with file IDs as labels
    ax.set_xticks(x_positions)
    ax.set_xticklabels(file_ids, rotation=45, ha='right')

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    unique_labels = {label: handle for label, handle in by_label.items()}
    legend = ax.legend(unique_labels.values(), unique_labels.keys(), title='Orientation', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Add sequence information below the legend for R1/R2 mode
    if use_separate_r1_r2:
        # Format sequences for display (truncate if too long)
        r1_display = [seq[:20] + '...' if len(seq) > 20 else seq for seq in r1_sequences]
        r2_display = [seq[:20] + '...' if len(seq) > 20 else seq for seq in r2_sequences]
        
        sequence_text = f"R1 sequences:\n{chr(10).join(r1_display)}\n\nR2 sequences:\n{chr(10).join(r2_display)}"
        
        # Position text below legend using axes coordinates
        ax.text(1.05, 0.65, sequence_text, transform=ax.transAxes, fontsize=8, 
                verticalalignment='top', horizontalalignment='left',
                bbox=dict(boxstyle="round,pad=0.5", facecolor="lightgray", alpha=0.8))

    ax.set_xlabel('FASTQ file type')
    ax.set_ylabel('Percentage of reads')
    
    if use_separate_r1_r2:
        ax.set_title('Read orientation analysis: Forward (R1) vs Reverse (R2) sequences')
    else:
        ax.set_title(f'Reads containing target sequences')
    
    plt.tight_layout()
    plt.savefig(plot_png_file)
    logging.info(f"Plot saved to {plot_png_file}")
    plt.savefig(plot_pdf_file)
    logging.info(f"Plot saved to {plot_pdf_file}")
    
    # Log absolute paths for debugging
    logging.info(f"=== FILES CREATED ===")
    logging.info(f"Results TSV: {os.path.abspath(results_file)}")
    logging.info(f"Plot PNG: {os.path.abspath(plot_png_file)}")
    logging.info(f"Plot PDF: {os.path.abspath(plot_pdf_file)}")
    logging.info(f"Working directory: {os.getcwd()}")
    logging.info(f"Output prefix used: {output_file_prefix}")
    logging.info(f"=====================")

def parse_args():
    """Parse command line arguments when script is called directly"""
    parser = argparse.ArgumentParser(description="Create read reorientation barplots")
    parser.add_argument('--fastq-files', nargs='+', help='FASTQ files to analyze')
    parser.add_argument('--sequences', nargs='+', help='Sequences to search for (combined mode)')
    parser.add_argument('--r1-sequences', nargs='+', help='R1 sequences to search for (separate R1/R2 mode)')
    parser.add_argument('--r2-sequences', nargs='+', help='R2 sequences to search for (separate R1/R2 mode)')
    parser.add_argument('--output-prefix', default=None, 
                        help='Output file prefix (default: uses config output_dir/reorientation_plot or ./reorientation_plot)')
    parser.add_argument('--n-reads', type=int, default=100000, 
                        help='Number of reads to analyze per file (default: 100000)')
    parser.add_argument('--config', help='Config file path (if not provided, will auto-detect)')
    parser.add_argument('--sample', help='Sample name for auto-detection of FASTQ files')
    return parser.parse_args()

def get_fastq_files_from_config(config, sample_name=None):
    """Get FASTQ files from config and metadata"""
    if not pc or not pm:
        raise ImportError("Config parsing utilities not available")
    
    # Get metadata file
    metadata_file = config.get('metadata_file', config.get('METADATA_FILE'))
    if not metadata_file:
        # Auto-detect metadata file
        if os.path.exists("sample_metadata.yaml"):
            metadata_file = "sample_metadata.yaml"
        elif os.path.exists("sample_metadata.txt"):
            metadata_file = "sample_metadata.txt"
        else:
            raise FileNotFoundError("No metadata file found")
    
    # Get output directory (where combined files should be)
    output_dir = config.get('output_dir', config.get('OUTPUT_DIR', 
                           config.get('docker_output_dir', config.get('DOCKER_OUTPUT_DIR', ''))))
    
    if not output_dir:
        raise ValueError("Output directory not specified in config")
    
    # Parse metadata to get sample/library information
    libraries, sample_names, fastq_files = pm.process_metadata(metadata_file, "")  # Don't need fastq_dir for this
    
    fastq_files_to_analyze = []
    samples_to_process = []
    
    if sample_name:
        # Process only specified sample
        if sample_name not in sample_names:
            raise ValueError(f"Sample '{sample_name}' not found in metadata")
        sample_idx = sample_names.index(sample_name)
        samples_to_process = [(libraries[sample_idx], sample_names[sample_idx])]
    else:
        # Process all samples
        samples_to_process = list(zip(libraries, sample_names))
    
    # Look for combined and reoriented files for each sample
    for library, sample in samples_to_process:
        # Expected file patterns for this library
        combined_r1 = os.path.join(output_dir, f"{library}_combined_R1.fq.gz")
        combined_r2 = os.path.join(output_dir, f"{library}_combined_R2.fq.gz")
        reoriented_r1 = os.path.join(output_dir, f"{library}_combined_R1.reoriented.fq.gz")
        reoriented_r2 = os.path.join(output_dir, f"{library}_combined_R2.reoriented.fq.gz")
        
        # Check which files exist and add them
        files_for_sample = []
        if os.path.exists(combined_r1):
            files_for_sample.append(combined_r1)
        if os.path.exists(combined_r2):
            files_for_sample.append(combined_r2)
        if os.path.exists(reoriented_r1):
            files_for_sample.append(reoriented_r1)
        if os.path.exists(reoriented_r2):
            files_for_sample.append(reoriented_r2)
        
        if not files_for_sample:
            logging.warning(f"No combined/reoriented files found for library {library} in {output_dir}")
        else:
            fastq_files_to_analyze.extend(files_for_sample)
            logging.info(f"Found {len(files_for_sample)} files for library {library}: {[os.path.basename(f) for f in files_for_sample]}")
    
    if not fastq_files_to_analyze:
        raise FileNotFoundError("No combined/reoriented FASTQ files found in output directory")
    
    # Return files and a sample identifier
    if sample_name:
        return fastq_files_to_analyze, sample_name
    else:
        return fastq_files_to_analyze, "all_samples"

def process_samples_from_config(config, sample_name, sequences, output_prefix, n_reads, r1_sequences=None, r2_sequences=None):
    """Process samples individually, creating separate plots and TSV files for each"""
    if not pc or not pm:
        raise ImportError("Config parsing utilities not available")
    
    # Get metadata file
    metadata_file = config.get('metadata_file', config.get('METADATA_FILE'))
    if not metadata_file:
        # Auto-detect metadata file
        if os.path.exists("sample_metadata.yaml"):
            metadata_file = "sample_metadata.yaml"
        elif os.path.exists("sample_metadata.txt"):
            metadata_file = "sample_metadata.txt"
        else:
            raise FileNotFoundError("No metadata file found")
    
    # Get output directory (where combined files should be)
    output_dir = config.get('output_dir', config.get('OUTPUT_DIR', 
                           config.get('docker_output_dir', config.get('DOCKER_OUTPUT_DIR', ''))))
    
    if not output_dir:
        raise ValueError("Output directory not specified in config")
    
    # Parse metadata to get sample/library information
    libraries, sample_names, fastq_files = pm.process_metadata(metadata_file, "")
    
    samples_to_process = []
    if sample_name:
        # Process only specified sample
        if sample_name not in sample_names:
            raise ValueError(f"Sample '{sample_name}' not found in metadata")
        sample_idx = sample_names.index(sample_name)
        samples_to_process = [(libraries[sample_idx], sample_names[sample_idx])]
    else:
        # Process all samples
        samples_to_process = list(zip(libraries, sample_names))
    
    # Summary data for all samples
    summary_data = []
    summary_headers = ["Library", "Sample", "File_Type", "Found_Percentage", "Not_Found_Percentage", 
                      "Total_Reads", "Found_Reads", "Primary_Seq", "Primary_Hits", "Primary_Percentage"]
    
    # Process each sample individually
    for library, sample in samples_to_process:
        logging.info(f"Processing library {library} (sample: {sample})")
        
        # Expected file patterns for this library
        combined_r1 = os.path.join(output_dir, f"{library}_combined_R1.fq.gz")
        combined_r2 = os.path.join(output_dir, f"{library}_combined_R2.fq.gz")
        reoriented_r1 = os.path.join(output_dir, f"{library}_combined_R1.reoriented.fq.gz")
        reoriented_r2 = os.path.join(output_dir, f"{library}_combined_R2.reoriented.fq.gz")
        
        # Check which files exist
        files_for_sample = []
        if os.path.exists(combined_r1):
            files_for_sample.append(combined_r1)
        if os.path.exists(combined_r2):
            files_for_sample.append(combined_r2)
        if os.path.exists(reoriented_r1):
            files_for_sample.append(reoriented_r1)
        if os.path.exists(reoriented_r2):
            files_for_sample.append(reoriented_r2)
        
        if not files_for_sample:
            logging.warning(f"No combined/reoriented files found for library {library} in {output_dir}")
            continue
        
        logging.info(f"Found {len(files_for_sample)} files for library {library}: {[os.path.basename(f) for f in files_for_sample]}")
        
        # Create sample-specific output prefix with both library and sample name
        # Apply the same output directory logic as main function
        if not output_prefix.startswith('/'):
            # Relative prefix, prepend with output directory
            docker_paths = config.get('docker_paths', False)
            if docker_paths:
                config_output_dir = config.get('docker_output_dir', config.get('output_dir', ''))
            else:
                config_output_dir = config.get('output_dir', config.get('docker_output_dir', ''))
            
            if config_output_dir:
                sample_output_prefix = os.path.join(config_output_dir, f"{output_prefix}_{library}_{sample}")
            else:
                sample_output_prefix = f"{output_prefix}_{library}_{sample}"
        else:
            # Absolute prefix, use as-is
            sample_output_prefix = f"{output_prefix}_{library}_{sample}"
        
        logging.info(f"Sample output prefix: {sample_output_prefix}")
        
        # Run analysis for this sample
        logging.info(f"Analyzing FASTQ files for {library} ({sample}): {files_for_sample}")
        plot_positions(files_for_sample, sequences, sample_output_prefix, n_reads, sample, library, r1_sequences, r2_sequences)
        
        # Extract summary data from the results file
        results_file = f"{sample_output_prefix}_results.tsv"
        if os.path.exists(results_file):
            sample_summary = extract_summary_from_results(results_file, library, sample)
            summary_data.extend(sample_summary)
    
    # Create overall summary report if processing multiple samples
    if len(samples_to_process) > 1:
        # Apply the same output directory logic for summary
        if not output_prefix.startswith('/'):
            # Relative prefix, prepend with output directory
            docker_paths = config.get('docker_paths', False)
            if docker_paths:
                config_output_dir = config.get('docker_output_dir', config.get('output_dir', ''))
            else:
                config_output_dir = config.get('output_dir', config.get('docker_output_dir', ''))
            
            if config_output_dir:
                summary_output_prefix = os.path.join(config_output_dir, output_prefix)
            else:
                summary_output_prefix = output_prefix
        else:
            # Absolute prefix, use as-is
            summary_output_prefix = output_prefix
        
        logging.info(f"Summary output prefix: {summary_output_prefix}")
        create_summary_report(summary_data, summary_headers, summary_output_prefix)

def extract_summary_from_results(results_file, library, sample):
    """Extract summary statistics from a results TSV file"""
    summary_rows = []
    
    try:
        with open(results_file, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            return summary_rows
            
        # Check if this is R1/R2 separate mode or combined mode
        header = lines[0].strip().split('\t')
        is_separate_mode = 'Sequence_Type' in header and 'R1_Hits' in header and 'R2_Hits' in header
        
        if is_separate_mode:
            # Handle R1/R2 separate mode
            file_data = {}
            for line in lines[1:]:  # Skip header
                parts = line.strip().split('\t')
                if len(parts) >= 13:  # New format has 13 columns
                    lib_col, sample_col, file_path, file_id, seq_type, sequence, position, count, seq_hits, total, r1_hits, r2_hits, any_seq_hits = parts[:13]
                    
                    if file_id not in file_data:
                        file_data[file_id] = {
                            'total_reads': int(total),
                            'r1_hits': int(r1_hits),
                            'r2_hits': int(r2_hits),
                            'any_seq_hits': int(any_seq_hits),
                            'r1_sequences': {},
                            'r2_sequences': {}
                        }
                    
                    # Track individual sequence hits
                    if seq_type == 'R1':
                        if sequence not in file_data[file_id]['r1_sequences']:
                            file_data[file_id]['r1_sequences'][sequence] = int(seq_hits)
                    elif seq_type == 'R2':
                        if sequence not in file_data[file_id]['r2_sequences']:
                            file_data[file_id]['r2_sequences'][sequence] = int(seq_hits)
            
            # Create summary rows for each file type
            for file_id, data in file_data.items():
                total_reads = data['total_reads']
                r1_hits = data['r1_hits']
                r2_hits = data['r2_hits']
                not_found = total_reads - r1_hits - r2_hits
                
                r1_percentage = (r1_hits / total_reads) * 100 if total_reads > 0 else 0
                r2_percentage = (r2_hits / total_reads) * 100 if total_reads > 0 else 0
                not_found_percentage = (not_found / total_reads) * 100 if total_reads > 0 else 0
                
                # Create base row
                base_row = [library, sample, file_id, f"{r1_percentage:.2f}", f"{r2_percentage:.2f}", 
                           f"{not_found_percentage:.2f}", total_reads, r1_hits, r2_hits, not_found]
                
                # Add individual R1 sequence columns
                r1_seq_cols = []
                for seq, hits in data['r1_sequences'].items():
                    percentage = (hits / total_reads) * 100 if total_reads > 0 else 0
                    r1_seq_cols.extend([seq, hits, f"{percentage:.2f}"])
                
                # Add individual R2 sequence columns  
                r2_seq_cols = []
                for seq, hits in data['r2_sequences'].items():
                    percentage = (hits / total_reads) * 100 if total_reads > 0 else 0
                    r2_seq_cols.extend([seq, hits, f"{percentage:.2f}"])
                
                summary_rows.append(base_row + r1_seq_cols + r2_seq_cols)
        
        else:
            # Handle combined mode (backward compatibility)
            file_data = {}
            for line in lines[1:]:  # Skip header
                parts = line.strip().split('\t')
                if len(parts) >= 10:  # Original format
                    library_col, sample_col, file_path, file_id, sequence, position, count, seq_hits, total, any_seq_hits = parts[:10]
                    
                    if file_id not in file_data:
                        file_data[file_id] = {
                            'total_reads': int(total),
                            'found_reads': int(any_seq_hits),
                            'sequences': {}
                        }
                    
                    if sequence not in file_data[file_id]['sequences']:
                        file_data[file_id]['sequences'][sequence] = int(seq_hits)
            
            # Create summary rows for each file type
            for file_id, data in file_data.items():
                total_reads = data['total_reads']
                found_reads = data['found_reads']
                found_percentage = (found_reads / total_reads) * 100 if total_reads > 0 else 0
                not_found_percentage = 100 - found_percentage
                
                # Find primary sequence (one with most hits)
                primary_seq = max(data['sequences'].items(), key=lambda x: x[1]) if data['sequences'] else ("", 0)
                primary_seq_name, primary_hits = primary_seq
                primary_percentage = (primary_hits / total_reads) * 100 if total_reads > 0 else 0
                
                summary_rows.append([
                    library, sample, file_id, f"{found_percentage:.2f}", f"{not_found_percentage:.2f}",
                    total_reads, found_reads, primary_seq_name, primary_hits, f"{primary_percentage:.2f}"
                ])
    
    except Exception as e:
        logging.error(f"Error extracting summary from {results_file}: {e}")
    
    return summary_rows

def create_summary_report(summary_data, headers, output_prefix):
    """Create an overall summary report for all samples"""
    summary_file = f"{output_prefix}_all_samples_summary.tsv"
    
    try:
        if not summary_data:
            logging.warning("No summary data to write")
            return
            
        # Determine if we're using R1/R2 separate mode based on the first row
        first_row = summary_data[0]
        is_separate_mode = len(first_row) > 10  # Separate mode has more columns
        
        if is_separate_mode:
            # Dynamic headers for R1/R2 mode - determine from the data
            base_headers = ["Library", "Sample", "File_Type", "R1_Found_Percentage", "R2_Found_Percentage", 
                          "Not_Found_Percentage", "Total_Reads", "R1_Hits", "R2_Hits", "Not_Found"]
            
            # Find all unique R1 and R2 sequences across all samples
            r1_sequences = set()
            r2_sequences = set()
            
            for row in summary_data:
                # Extract sequence information from the variable part of each row
                variable_part = row[10:]  # Everything after the base columns
                
                # Process in groups of 3 (sequence, hits, percentage)
                for i in range(0, len(variable_part), 3):
                    if i + 2 < len(variable_part):
                        seq_name = variable_part[i]
                        # Determine if it's R1 or R2 based on context or naming
                        # We'll need to get this info from the original analysis
                        pass
            
            # For now, use a simpler approach - create dynamic headers based on actual data
            max_cols = max(len(row) for row in summary_data)
            dynamic_headers = base_headers.copy()
            
            # Add placeholder headers for the additional sequence columns
            extra_cols = max_cols - len(base_headers)
            for i in range(0, extra_cols, 3):
                seq_num = (i // 3) + 1
                dynamic_headers.extend([f"Sequence_{seq_num}_Name", f"Sequence_{seq_num}_Hits", f"Sequence_{seq_num}_Percentage"])
            
            actual_headers = dynamic_headers[:max_cols]
        else:
            # Use provided headers for combined mode
            actual_headers = headers
        
        with open(summary_file, 'w') as f:
            f.write('\t'.join(actual_headers) + '\n')
            for row in summary_data:
                # Pad row if necessary to match header length
                padded_row = row + [''] * (len(actual_headers) - len(row))
                f.write('\t'.join(map(str, padded_row[:len(actual_headers)])) + '\n')
        
        logging.info(f"Summary report saved to {summary_file}")
        
        # Log absolute path for debugging
        logging.info(f"=== SUMMARY FILE CREATED ===")
        logging.info(f"Summary TSV: {os.path.abspath(summary_file)}")
        logging.info(f"Working directory: {os.getcwd()}")
        logging.info(f"Output prefix used: {output_prefix}")
        logging.info(f"============================")
        
        # Only show console summary if we have meaningful data
        logging.info("=== REORIENTATION SUMMARY ===")
        logging.info(f"Processed {len(set([row[0] for row in summary_data]))} libraries")
        logging.info(f"Detailed results available in {summary_file}")
        
        if is_separate_mode:
            logging.info("Summary includes R1/R2 separate analysis with individual sequence statistics")
        else:
            logging.info("Summary includes combined sequence analysis")
    
    except Exception as e:
        logging.error(f"Error creating summary report: {e}")

def main():
    """Main function when script is called directly"""
    args = parse_args()
    
    # Initialize config to None
    config = None
    
    # If fastq files and sequences are provided directly, try to load config for output directory
    if args.fastq_files and (args.sequences or (args.r1_sequences and args.r2_sequences)):
        logging.info("Using provided FASTQ files and sequences")
        
        # Try to load config for output directory even when using direct parameters
        if pc:
            config_file = args.config
            if not config_file:
                if os.path.exists("config.yaml"):
                    config_file = "config.yaml"
                elif os.path.exists("config.sh"):
                    config_file = "config.sh"
            
            if config_file:
                try:
                    logging.info(f"Loading config from {config_file} for output directory")
                    config = pc.parse_config(config_file, quiet=True)
                except Exception as e:
                    logging.warning(f"Could not load config file {config_file}: {e}")
        
        # Determine output prefix - use config output directory if available
        if args.output_prefix is None:
            # No output prefix specified, use config output directory with default name
            if config:
                # Check if using docker paths
                docker_paths = config.get('docker_paths', False)
                if docker_paths:
                    output_dir = config.get('docker_output_dir', config.get('output_dir', ''))
                else:
                    output_dir = config.get('output_dir', config.get('docker_output_dir', ''))
                
                if output_dir:
                    output_prefix = os.path.join(output_dir, 'reorientation_plot')
                    logging.info(f"Using default output prefix in config output directory: {output_prefix}")
                else:
                    output_prefix = 'reorientation_plot'
                    logging.info(f"No output dir in config, using default: {output_prefix}")
            else:
                output_prefix = 'reorientation_plot'
                logging.info(f"No config available, using default: {output_prefix}")
        elif not args.output_prefix.startswith('/'):
            # Relative output prefix specified, use config output directory if available
            if config:
                # Check if using docker paths
                docker_paths = config.get('docker_paths', False)
                if docker_paths:
                    output_dir = config.get('docker_output_dir', config.get('output_dir', ''))
                else:
                    output_dir = config.get('output_dir', config.get('docker_output_dir', ''))
                
                if output_dir:
                    output_prefix = os.path.join(output_dir, args.output_prefix)
                    logging.info(f"Using output directory from config: {output_dir}")
                else:
                    output_prefix = args.output_prefix
                    logging.info(f"No output dir in config, using prefix as-is: {output_prefix}")
            else:
                output_prefix = args.output_prefix
                logging.info(f"No config available, using prefix as-is: {output_prefix}")
        else:
            # Absolute path specified, use as-is
            output_prefix = args.output_prefix
        
        if args.r1_sequences and args.r2_sequences:
            plot_positions(args.fastq_files, None, output_prefix, args.n_reads, None, None, args.r1_sequences, args.r2_sequences)
        else:
            plot_positions(args.fastq_files, args.sequences, output_prefix, args.n_reads, None, None, None, None)
        return
    
    # Otherwise, parse config file (always try to parse config for missing parameters)
    if not pc:
        if args.fastq_files and (args.sequences or (args.r1_sequences and args.r2_sequences)):
            # We have what we need, proceed without config (output to current directory)
            logging.warning("Config parsing utilities not available - using current directory for output")
            output_prefix = args.output_prefix if args.output_prefix is not None else 'reorientation_plot'
            if args.r1_sequences and args.r2_sequences:
                plot_positions(args.fastq_files, None, output_prefix, args.n_reads, None, None, args.r1_sequences, args.r2_sequences)
            else:
                plot_positions(args.fastq_files, args.sequences, output_prefix, args.n_reads, None, None, None, None)
            return
        else:
            raise ImportError("Config parsing utilities not available and missing required parameters. Please provide --fastq-files and --sequences (or --r1-sequences and --r2-sequences) directly.")
    
    # Auto-detect or use provided config file (matches other pipeline scripts)
    config_file = args.config
    if not config_file:
        if os.path.exists("config.yaml"):
            config_file = "config.yaml"
        elif os.path.exists("config.sh"):
            config_file = "config.sh"
        else:
            if args.fastq_files and (args.sequences or (args.r1_sequences and args.r2_sequences)):
                # We have explicit parameters, proceed without config (output to current directory)
                logging.info("No config file found, using provided parameters and current directory for output")
                output_prefix = args.output_prefix if args.output_prefix is not None else 'reorientation_plot'
                if args.r1_sequences and args.r2_sequences:
                    plot_positions(args.fastq_files, None, output_prefix, args.n_reads, None, None, args.r1_sequences, args.r2_sequences)
                else:
                    plot_positions(args.fastq_files, args.sequences, output_prefix, args.n_reads, None, None, None, None)
                return
            else:
                raise FileNotFoundError("No config file found and missing required parameters. Please specify --config or provide --fastq-files and --sequences (or --r1-sequences and --r2-sequences)")
    
    logging.info(f"Loading config from {config_file}")
    config = pc.parse_config(config_file, quiet=True)
    
    # Get search sequences from config or command line
    if args.r1_sequences and args.r2_sequences:
        # Use command line R1/R2 sequences
        sequences = None
        r1_sequences = args.r1_sequences
        r2_sequences = args.r2_sequences
    elif args.sequences:
        sequences = args.sequences
        r1_sequences = None
        r2_sequences = None
    else:
        # Get search sequences from config
        r1_seqs = config.get('r1_seqs', config.get('R1_SEQS', []))
        r2_seqs = config.get('r2_seqs', config.get('R2_SEQS', []))
        
        if isinstance(r1_seqs, str):
            r1_seqs = [r1_seqs]
        if isinstance(r2_seqs, str):
            r2_seqs = [r2_seqs]
        
        # Check if we have both R1 and R2 sequences for separate analysis
        if r1_seqs and r2_seqs:
            # Use separate R1/R2 mode
            r1_sequences = r1_seqs
            r2_sequences = r2_seqs
            sequences = None
            logging.info("Using separate R1/R2 sequence analysis mode")
        else:
            # Use combined mode (backward compatibility)
            if r1_seqs:
                sequences = r1_seqs.copy()
                if r2_seqs:
                    logging.warning("R2 sequences found but no R1 sequences - ignoring R2 sequences for backward compatibility")
            elif r2_seqs:
                sequences = r2_seqs
            else:
                raise ValueError("No search sequences found in config file. Please provide --sequences or ensure r1_seqs/R1_SEQS is defined in config")
            r1_sequences = None
            r2_sequences = None
            logging.info("Using combined sequence analysis mode (R2 sequences ignored)")
    
    # Get FASTQ files from command line or config
    if args.fastq_files:
        fastq_files = args.fastq_files
        
        # Determine output prefix - use config output directory if available
        if args.output_prefix is None:
            # No output prefix specified, use config output directory with default name
            if config:
                # Check if using docker paths
                docker_paths = config.get('docker_paths', False)
                if docker_paths:
                    output_dir = config.get('docker_output_dir', config.get('output_dir', ''))
                else:
                    output_dir = config.get('output_dir', config.get('docker_output_dir', ''))
                
                if output_dir:
                    output_prefix = os.path.join(output_dir, 'reorientation_plot')
                else:
                    output_prefix = 'reorientation_plot'
                logging.info(f"Calculated output prefix (fastq files branch): {output_prefix}")
            else:
                output_prefix = 'reorientation_plot'
                logging.info(f"No config available (fastq files branch): {output_prefix}")
        elif not args.output_prefix.startswith('/'):
            # Relative output prefix specified, use config output directory if available
            if config:
                # Check if using docker paths
                docker_paths = config.get('docker_paths', False)
                if docker_paths:
                    output_dir = config.get('docker_output_dir', config.get('output_dir', ''))
                else:
                    output_dir = config.get('output_dir', config.get('docker_output_dir', ''))
                
                if output_dir:
                    output_prefix = os.path.join(output_dir, args.output_prefix)
                else:
                    output_prefix = args.output_prefix
                logging.info(f"Calculated output prefix (fastq files branch): {output_prefix}")
            else:
                output_prefix = args.output_prefix
                logging.info(f"No config available (fastq files branch): {output_prefix}")
        else:
            # Absolute path specified, use as-is
            output_prefix = args.output_prefix
            
        # Run single analysis for provided files
        logging.info(f"Analyzing FASTQ files: {fastq_files}")
        if sequences:
            logging.info(f"Searching for sequences: {sequences}")
        else:
            logging.info(f"Searching for R1 sequences: {r1_sequences} and R2 sequences: {r2_sequences}")
        plot_positions(fastq_files, sequences, output_prefix, args.n_reads, None, None, r1_sequences, r2_sequences)
    else:
        try:
            # Process samples individually when using config
            # Use default prefix if none specified
            output_prefix = args.output_prefix if args.output_prefix is not None else 'reorientation_plot'
            process_samples_from_config(config, args.sample, sequences, output_prefix, args.n_reads, r1_sequences, r2_sequences)
        except Exception as e:
            raise ValueError(f"Could not determine FASTQ files from config: {e}. Please provide --fastq-files explicitly")

if __name__ == "__main__":
    main()
