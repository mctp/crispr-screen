#!/usr/bin/env python

"""
view_fastq.py

- Views FASTQ file sequences with optional quality-based masking
- Can process individual files or use config/metadata for automatic file selection
- Supports quality threshold masking to hide low-quality bases

Dependencies:
    pip install biopython pyyaml pandas
"""

from html import parser
import sys
import os
import argparse
import gzip
from Bio import SeqIO
from pathlib import Path
import logging

# project-specific imports
import utilities.process_metadata as pm
import utilities.parse_config as pc

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    return logging.getLogger(__name__)

def parse_args():
    parser = argparse.ArgumentParser(description="View FASTQ file sequences with optional quality masking")
    parser.add_argument('-i', '--in-fastq', '--in-fq', '--fastq', '--fq', 
                        help='Input FASTQ file', type=Path)
    parser.add_argument('-s', '--sample', help='Sample name (if not provided, uses first sample)')
    parser.add_argument('--config', help='Config file path (config.yaml or config.sh)', default='config.sh')
    parser.add_argument('--metadata', help='Metadata file path (sample_metadata.yaml or sample_metadata.txt)', 
                        default='sample_metadata.txt')
    parser.add_argument('-n', '--nreads', type=int, default=10, 
                        help='Number of reads to display (default: 10)')
    parser.add_argument('-r', '--random', '--sample-random', action='store_true',
                        help='Randomly sample reads instead of taking first n reads (requires seqtk)')
    parser.add_argument('-q', '--quality-threshold', nargs='?', const=30, type=int, default=None,
                        help='Quality threshold for masking bases (Phred score). Use -q alone for default (30), -q VALUE for custom, or omit for no masking')
    parser.add_argument('-Q', '--avg-quality-threshold', nargs='?', const=25.0, type=float, default=None,
                        help='Average quality threshold for filtering reads (Phred score). Use -Q alone for default (25.0), -Q VALUE for custom, or omit for no filtering')
    parser.add_argument('--start-read', type=int, default=1,
                        help='Starting read number (1-based, default: 1)')
    parser.add_argument('--read-id', type=str, default=None,
                        help='Display specific read by ID')
    parser.add_argument('--show-quality', action='store_true',
                        help='Show quality scores alongside sequences')
    parser.add_argument('--show-header', action='store_true',
                        help='Show FASTQ headers')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Show detailed information (default: only show sequences)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for sampling (default: 42)')
    return parser.parse_args()

def get_parameters(logger):
    """Get parameters from args, config, or metadata parsing"""
    args = parse_args()
    
    # Check if required args are provided directly
    if args.in_fastq:
        return args.in_fastq, args
    
    # Load config file
    logger.info(f"Loading config from {args.config}")
    try:
        config = pc.parse_config(args.config, quiet=True)
    except Exception as e:
        logger.error(f"Failed to parse config file {args.config}: {e}")
        sys.exit(1)
    
    # Get paths from config
    metadata_file = args.metadata or config.get('metadata_file', config.get('METADATA_FILE', 'sample_metadata.txt'))
    fastq_dir = config.get('fastq_dir', config.get('FASTQ_DIR', config.get('docker_fastq_dir', config.get('DOCKER_FASTQ_DIR', ''))))
    
    logger.info(f"Using metadata file: {metadata_file}")
    logger.info(f"Using fastq directory: {fastq_dir}")
    
    # Parse metadata
    logger.info("Parsing metadata...")
    try:
        libraries, sample_names, fastq_files = pm.process_metadata(metadata_file, fastq_dir)
        
        if not libraries:
            raise ValueError("No samples found in metadata file")
        
        # Use specified sample or first sample
        if args.sample:
            # Find the specified sample
            try:
                sample_idx = sample_names.index(args.sample)
                sample = args.sample
                fastq_file_list = fastq_files[sample_idx]
            except ValueError:
                raise ValueError(f"Sample '{args.sample}' not found in metadata")
        else:
            # Use first sample
            sample_idx = 0
            sample = sample_names[sample_idx]
            fastq_file_list = fastq_files[sample_idx]
            logger.info(f"Using first sample from metadata: {sample}")
        
        # Parse fastq file path - use R1 file
        if ';' in fastq_file_list:
            r1_files, r2_files = fastq_file_list.split(';')
            fastq_file = Path(r1_files.split(',')[0])  # Use first R1 file
        else:
            fastq_file = Path(fastq_file_list.split(',')[0])  # Use first file
        
        logger.info(f"Using FASTQ file: {fastq_file}")
        logger.info(f"Using sample: {sample}")
        
    except Exception as e:
        logger.error(f"Failed to parse metadata: {e}")
        sys.exit(1)
    
    return fastq_file, args

def check_tool(tool):
    """Check if required tool is available"""
    from shutil import which
    if not which(tool):
        return False
    return True

def sample_fastq_with_seqtk(fastq_file, n_reads, seed=42):
    """Use seqtk to randomly sample reads from FASTQ file"""
    import tempfile
    import subprocess
    
    # Create temporary file for sampled reads
    temp_file = tempfile.NamedTemporaryFile(mode='w+', suffix='.fastq', delete=False)
    temp_path = temp_file.name
    temp_file.close()
    
    try:
        # Build seqtk command
        if str(fastq_file).endswith('.gz'):
            cmd = f"seqtk sample -s {seed} {fastq_file} {n_reads} > {temp_path}"
        else:
            cmd = f"seqtk sample -s {seed} {fastq_file} {n_reads} > {temp_path}"
        
        # Run seqtk
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            raise RuntimeError(f"seqtk failed: {result.stderr}")
        
        return temp_path
        
    except Exception as e:
        # Clean up temp file on error
        if os.path.exists(temp_path):
            os.unlink(temp_path)
        raise e

def open_fastq(file_path):
    """Open FASTQ file, handling both gzipped and plain text files"""
    if str(file_path).endswith('.gz'):
        return gzip.open(file_path, 'rt')
    else:
        return open(file_path, 'r')

def mask_low_quality_bases(sequence, quality_scores, threshold):
    """
    Mask bases below quality threshold with colored background while preserving spacing
    """
    if threshold is None:
        return sequence
    
    masked_sequence = ""
    for base, qual in zip(sequence, quality_scores):
        # Convert quality character to Phred score (assuming Phred+33 encoding)
        phred_score = ord(qual) - 33
        if phred_score < threshold:
            masked_sequence += f"\033[93m\033[41m{base}\033[0m"  # Bright yellow text on red background
        else:
            masked_sequence += base
    
    return masked_sequence

def format_quality_scores(quality_scores, threshold=None):
    """Format quality scores for display, highlighting low quality"""
    if threshold is None:
        return quality_scores
    
    formatted = ""
    for qual in quality_scores:
        phred_score = ord(qual) - 33
        if phred_score < threshold:
            formatted += f"\033[91m{qual}\033[0m"  # Red for low quality
        else:
            formatted += qual
    
    return formatted

    """View FASTQ sequences with optional quality masking"""
    logger = logging.getLogger(__name__)
    
    if not fastq_file.exists():
        if args.verbose:
            logger.error(f"FASTQ file does not exist: {fastq_file}")
        sys.exit(1)
    
    if args.verbose:
        logger.info(f"Reading FASTQ file: {fastq_file}")
        if args.quality_threshold is not None:
            logger.info(f"Using quality threshold: {args.quality_threshold} (Phred score)")
        else:
            logger.info("No quality masking (use -q for default threshold of 30)")
        logger.info(f"Displaying {args.nreads} reads starting from read {args.start_read}")
    
    read_count = 0
    displayed_count = 0
    
    try:
        with open_fastq(fastq_file) as handle:
            for record in SeqIO.parse(handle, "fastq"):
                read_count += 1
                
                # Skip reads until we reach the starting read
                if read_count < args.start_read:
                    continue
                
                # Check if we're looking for a specific read ID
                if args.read_id and args.read_id not in record.id:
                    continue
                
                # Stop if we've displayed enough reads
                if not args.read_id and displayed_count >= args.nreads:
                    break
                
                displayed_count += 1
                
                # Get sequence and quality
                sequence = str(record.seq)
                quality_scores = record.letter_annotations.get("phred_quality", None)
                
                if quality_scores is None:
                    # Use raw quality string if phred_quality not available
                    quality_string = record.format("fastq").split('\n')[3]
                else:
                    # Convert phred scores back to quality string
                    quality_string = ''.join(chr(q + 33) for q in quality_scores)
                
                # Apply quality masking if threshold is set
                if args.quality_threshold is not None:
                    masked_sequence = mask_low_quality_bases(sequence, quality_string, args.quality_threshold)
                    display_sequence = masked_sequence
                else:
                    display_sequence = sequence
                
                # VERBOSE MODE: Show detailed information
                if args.verbose:
                    print(f"\n{'='*80}")
                    print(f"Read {read_count}: {record.id}")
                    
                    if args.show_header:
                        print(f"Header: {record.description}")
                    
                    if args.quality_threshold is not None:
                        print(f"Sequence: {display_sequence}")
                        print(f"Original: {sequence}")
                    else:
                        print(f"Sequence: {display_sequence}")
                    
                    print(f"Length:   {len(sequence)} bp")
                    
                    # Show quality scores if requested
                    if args.show_quality:
                        if args.quality_threshold is not None:
                            formatted_quality = format_quality_scores(quality_string, args.quality_threshold)
                            print(f"Quality:  {formatted_quality}")
                        else:
                            print(f"Quality:  {quality_string}")
                        
                        # Show quality statistics
                        if quality_scores is not None:
                            avg_quality = sum(quality_scores) / len(quality_scores)
                            min_quality = min(quality_scores)
                            max_quality = max(quality_scores)
                            print(f"Quality stats: avg={avg_quality:.1f}, min={min_quality}, max={max_quality}")
                        else:
                            # Calculate from quality string
                            phred_scores = [ord(q) - 33 for q in quality_string]
                            avg_quality = sum(phred_scores) / len(phred_scores)
                            min_quality = min(phred_scores)
                            max_quality = max(phred_scores)
                            print(f"Quality stats: avg={avg_quality:.1f}, min={min_quality}, max={max_quality}")
                    
                    # Show masking statistics if threshold is used
                    if args.quality_threshold is not None:
                        phred_scores = [ord(q) - 33 for q in quality_string]
                        low_quality_count = sum(1 for q in phred_scores if q < args.quality_threshold)
                        masking_percentage = (low_quality_count / len(sequence)) * 100
                        print(f"Masked:   {low_quality_count}/{len(sequence)} bases ({masking_percentage:.1f}%)")
                
                # QUIET MODE: Only show the sequence
                else:
                    print(display_sequence)
                
                # If looking for specific read ID, stop after finding it
                if args.read_id and args.read_id in record.id:
                    break
    
    except Exception as e:
        if args.verbose:
            logger.error(f"Error reading FASTQ file: {e}")
        sys.exit(1)
    
    if args.verbose:
        print(f"\n{'='*80}")
        print(f"Displayed {displayed_count} reads from {fastq_file}")
        if args.quality_threshold is not None:
            print(f"Used quality threshold: {args.quality_threshold} (bases below this are masked)")
    """View FASTQ sequences with optional quality masking"""
    logger = logging.getLogger(__name__)
    
    if not fastq_file.exists():
        logger.error(f"FASTQ file does not exist: {fastq_file}")
        sys.exit(1)
    
    logger.info(f"Reading FASTQ file: {fastq_file}")
    if args.quality_threshold is not None:
        logger.info(f"Using quality threshold: {args.quality_threshold} (Phred score)")
    else:
        logger.info("No quality masking (use -q for default threshold of 30)")
    logger.info(f"Displaying {args.nreads} reads starting from read {args.start_read}")
    
    read_count = 0
    displayed_count = 0
    
    try:
        with open_fastq(fastq_file) as handle:
            for record in SeqIO.parse(handle, "fastq"):
                read_count += 1
                
                # Skip reads until we reach the starting read
                if read_count < args.start_read:
                    continue
                
                # Check if we're looking for a specific read ID
                if args.read_id and args.read_id not in record.id:
                    continue
                
                # Stop if we've displayed enough reads
                if not args.read_id and displayed_count >= args.nreads:
                    break
                
                displayed_count += 1
                
                # Display read information
                print(f"\n{'='*80}")
                print(f"Read {read_count}: {record.id}")
                
                if args.show_header:
                    print(f"Header: {record.description}")
                
                # Get sequence and quality
                sequence = str(record.seq)
                quality_scores = record.letter_annotations.get("phred_quality", None)
                
                if quality_scores is None:
                    # Use raw quality string if phred_quality not available
                    quality_string = record.format("fastq").split('\n')[3]
                else:
                    # Convert phred scores back to quality string
                    quality_string = ''.join(chr(q + 33) for q in quality_scores)
                
                # Apply quality masking if threshold is set
                if args.quality_threshold is not None:
                    masked_sequence = mask_low_quality_bases(sequence, quality_string, args.quality_threshold)
                    print(f"Sequence: {masked_sequence}")
                    print(f"Original: {sequence}")
                else:
                    print(f"Sequence: {sequence}")
                
                print(f"Length:   {len(sequence)} bp")
                
                # Show quality scores if requested
                if args.show_quality:
                    if args.quality_threshold is not None:
                        formatted_quality = format_quality_scores(quality_string, args.quality_threshold)
                        print(f"Quality:  {formatted_quality}")
                    else:
                        print(f"Quality:  {quality_string}")
                    
                    # Show quality statistics
                    if quality_scores is not None:
                        avg_quality = sum(quality_scores) / len(quality_scores)
                        min_quality = min(quality_scores)
                        max_quality = max(quality_scores)
                        print(f"Quality stats: avg={avg_quality:.1f}, min={min_quality}, max={max_quality}")
                    else:
                        # Calculate from quality string
                        phred_scores = [ord(q) - 33 for q in quality_string]
                        avg_quality = sum(phred_scores) / len(phred_scores)
                        min_quality = min(phred_scores)
                        max_quality = max(phred_scores)
                        print(f"Quality stats: avg={avg_quality:.1f}, min={min_quality}, max={max_quality}")
                
                # Show masking statistics if threshold is used
                if args.quality_threshold is not None:
                    phred_scores = [ord(q) - 33 for q in quality_string]
                    low_quality_count = sum(1 for q in phred_scores if q < args.quality_threshold)
                    masking_percentage = (low_quality_count / len(sequence)) * 100
                    print(f"Masked:   {low_quality_count}/{len(sequence)} bases ({masking_percentage:.1f}%)")
                
                # If looking for specific read ID, stop after finding it
                if args.read_id and args.read_id in record.id:
                    break
    
    except Exception as e:
        logger.error(f"Error reading FASTQ file: {e}")
        sys.exit(1)
    
    print(f"\n{'='*80}")
    print(f"Displayed {displayed_count} reads from {fastq_file}")
    if args.quality_threshold is not None:
        print(f"Used quality threshold: {args.quality_threshold} (bases below this are masked with ·)")

    """View FASTQ sequences with optional quality masking"""
    logger = logging.getLogger(__name__)
    
    if not fastq_file.exists():
        if args.verbose:
            logger.error(f"FASTQ file does not exist: {fastq_file}")
        sys.exit(1)
    
    if args.verbose:
        logger.info(f"Reading FASTQ file: {fastq_file}")
        if args.quality_threshold is not None:
            logger.info(f"Using quality threshold: {args.quality_threshold} (Phred score)")
        else:
            logger.info("No quality masking (use -q for default threshold of 30)")
        logger.info(f"Displaying {args.nreads} reads starting from read {args.start_read}")
    
    read_count = 0
    displayed_count = 0
    
    try:
        with open_fastq(fastq_file) as handle:
            for record in SeqIO.parse(handle, "fastq"):
                read_count += 1
                
                # Skip reads until we reach the starting read
                if read_count < args.start_read:
                    continue
                
                # Check if we're looking for a specific read ID
                if args.read_id and args.read_id not in record.id:
                    continue
                
                # Stop if we've displayed enough reads
                if not args.read_id and displayed_count >= args.nreads:
                    break
                
                displayed_count += 1
                
                # Get sequence and quality
                sequence = str(record.seq)
                quality_scores = record.letter_annotations.get("phred_quality", None)
                
                if quality_scores is None:
                    # Use raw quality string if phred_quality not available
                    quality_string = record.format("fastq").split('\n')[3]
                else:
                    # Convert phred scores back to quality string
                    quality_string = ''.join(chr(q + 33) for q in quality_scores)
                
                # Apply quality masking if threshold is set
                if args.quality_threshold is not None:
                    masked_sequence = mask_low_quality_bases(sequence, quality_string, args.quality_threshold)
                    display_sequence = masked_sequence
                else:
                    display_sequence = sequence
                
                # VERBOSE MODE: Show detailed information with formatting
                if args.verbose:
                    print(f"\n{'='*80}")
                    print(f"Read {read_count}: {record.id}")
                    
                    if args.show_header:
                        print(f"Header: {record.description}")
                    
                    if args.quality_threshold is not None:
                        print(f"Sequence: {display_sequence}")
                        print(f"Original: {sequence}")
                    else:
                        print(f"Sequence: {display_sequence}")
                    
                    print(f"Length:   {len(sequence)} bp")
                    
                    # Show quality scores if requested
                    if args.show_quality:
                        if args.quality_threshold is not None:
                            formatted_quality = format_quality_scores(quality_string, args.quality_threshold)
                            print(f"Quality:  {formatted_quality}")
                        else:
                            print(f"Quality:  {quality_string}")
                        
                        # Show quality statistics
                        if quality_scores is not None:
                            avg_quality = sum(quality_scores) / len(quality_scores)
                            min_quality = min(quality_scores)
                            max_quality = max(quality_scores)
                            print(f"Quality stats: avg={avg_quality:.1f}, min={min_quality}, max={max_quality}")
                        else:
                            # Calculate from quality string
                            phred_scores = [ord(q) - 33 for q in quality_string]
                            avg_quality = sum(phred_scores) / len(phred_scores)
                            min_quality = min(phred_scores)
                            max_quality = max(phred_scores)
                            print(f"Quality stats: avg={avg_quality:.1f}, min={min_quality}, max={max_quality}")
                    
                    # Show masking statistics if threshold is used
                    if args.quality_threshold is not None:
                        phred_scores = [ord(q) - 33 for q in quality_string]
                        low_quality_count = sum(1 for q in phred_scores if q < args.quality_threshold)
                        masking_percentage = (low_quality_count / len(sequence)) * 100
                        print(f"Masked:   {low_quality_count}/{len(sequence)} bases ({masking_percentage:.1f}%)")
                
                # QUIET MODE: Only show sequence (and quality if requested), no formatting
                else:
                    print(display_sequence)
                    if args.show_quality:
                        if args.quality_threshold is not None:
                            formatted_quality = format_quality_scores(quality_string, args.quality_threshold)
                            print(formatted_quality)
                        else:
                            print(quality_string)
                
                # If looking for specific read ID, stop after finding it
                if args.read_id and args.read_id in record.id:
                    break
    
    except Exception as e:
        if args.verbose:
            logger.error(f"Error reading FASTQ file: {e}")
        sys.exit(1)
    
    if args.verbose:
        print(f"\n{'='*80}")
        print(f"Displayed {displayed_count} reads from {fastq_file}")
        if args.quality_threshold is not None:
            print(f"Used quality threshold: {args.quality_threshold} (bases below this are masked)")

def calculate_average_quality(quality_string):
    """Calculate average Phred quality score from quality string"""
    if not quality_string:
        return 0
    phred_scores = [ord(q) - 33 for q in quality_string]
    return sum(phred_scores) / len(phred_scores)

def view_fastq_sequences(fastq_file, args):
    """View FASTQ sequences with optional quality masking"""
    logger = logging.getLogger(__name__)
    
    if not fastq_file.exists():
        if args.verbose:
            logger.error(f"FASTQ file does not exist: {fastq_file}")
        sys.exit(1)
    
    # Handle random sampling with adaptive quality filtering
    temp_file = None
    original_fastq = fastq_file  # Keep reference to original file
    
    if args.random:
        if not check_tool('seqtk'):
            if args.verbose:
                logger.error("seqtk not found in PATH. Random sampling requires seqtk.")
            else:
                print("ERROR: seqtk not found in PATH", file=sys.stderr)
            sys.exit(1)
        
        # If quality filtering is enabled, use adaptive sampling
        if args.avg_quality_threshold is not None:
            temp_file = adaptive_sample_with_quality_filter(original_fastq, args)
            fastq_file = Path(temp_file)
        else:
            # Simple random sampling without quality filtering
            if args.verbose:
                logger.info(f"Randomly sampling {args.nreads} reads using seqtk (seed: {args.seed})")
            
            try:
                temp_file = sample_fastq_with_seqtk(original_fastq, args.nreads, args.seed)
                fastq_file = Path(temp_file)
            except Exception as e:
                if args.verbose:
                    logger.error(f"Failed to sample reads: {e}")
                else:
                    print(f"ERROR: Failed to sample reads: {e}", file=sys.stderr)
                sys.exit(1)
    
    if args.verbose:
        logger.info(f"Reading FASTQ file: {fastq_file}")
        if args.quality_threshold is not None:
            logger.info(f"Using quality threshold for masking: {args.quality_threshold} (Phred score)")
        else:
            logger.info("No quality masking (use -q for default threshold of 30)")
        
        if args.avg_quality_threshold is not None:
            logger.info(f"Using average quality threshold for filtering: {args.avg_quality_threshold} (Phred score)")
        else:
            logger.info("No average quality filtering (use -Q for default threshold of 25)")
        
        if args.random:
            logger.info(f"Target: {args.nreads} reads after filtering")
        else:
            logger.info(f"Displaying {args.nreads} reads starting from read {args.start_read}")
    
    read_count = 0
    displayed_count = 0
    filtered_count = 0  # Track reads filtered by average quality
    target_reads = args.nreads
    
    try:
        with open_fastq(fastq_file) as handle:
            for record in SeqIO.parse(handle, "fastq"):
                read_count += 1
                
                # Skip reads until we reach the starting read (only if not random sampling)
                if not args.random and read_count < args.start_read:
                    continue
                
                # Check if we're looking for a specific read ID
                if args.read_id and args.read_id not in record.id:
                    continue
                
                # Get sequence and quality
                sequence = str(record.seq)
                quality_scores = record.letter_annotations.get("phred_quality", None)
                
                if quality_scores is None:
                    # Use raw quality string if phred_quality not available
                    quality_string = record.format("fastq").split('\n')[3]
                else:
                    # Convert phred scores back to quality string
                    quality_string = ''.join(chr(q + 33) for q in quality_scores)
                
                # Apply average quality filtering if threshold is set (and not using adaptive sampling)
                if args.avg_quality_threshold is not None and not args.random:
                    avg_quality = calculate_average_quality(quality_string)
                    if avg_quality < args.avg_quality_threshold:
                        filtered_count += 1
                        continue  # Skip this read
                
                # Stop if we've displayed enough reads
                if not args.read_id and displayed_count >= target_reads:
                    break
                
                displayed_count += 1
                
                # Apply quality masking if threshold is set
                if args.quality_threshold is not None:
                    masked_sequence = mask_low_quality_bases(sequence, quality_string, args.quality_threshold)
                    display_sequence = masked_sequence
                else:
                    display_sequence = sequence
                
                # VERBOSE MODE: Show detailed information with formatting
                if args.verbose:
                    print(f"\n{'='*80}")
                    if args.random:
                        print(f"Random Read {displayed_count}: {record.id}")
                    else:
                        print(f"Read {read_count}: {record.id}")
                    
                    if args.show_header:
                        print(f"Header: {record.description}")
                    
                    if args.quality_threshold is not None:
                        print(f"Sequence: {display_sequence}")
                        print(f"Original: {sequence}")
                    else:
                        print(f"Sequence: {display_sequence}")
                    
                    print(f"Length:   {len(sequence)} bp")
                    
                    # Show quality scores if requested
                    if args.show_quality:
                        if args.quality_threshold is not None:
                            formatted_quality = format_quality_scores(quality_string, args.quality_threshold)
                            print(f"Quality:  {formatted_quality}")
                        else:
                            print(f"Quality:  {quality_string}")
                        
                        # Show quality statistics
                        if quality_scores is not None:
                            avg_quality = sum(quality_scores) / len(quality_scores)
                            min_quality = min(quality_scores)
                            max_quality = max(quality_scores)
                            print(f"Quality stats: avg={avg_quality:.1f}, min={min_quality}, max={max_quality}")
                        else:
                            # Calculate from quality string
                            phred_scores = [ord(q) - 33 for q in quality_string]
                            avg_quality = sum(phred_scores) / len(phred_scores)
                            min_quality = min(phred_scores)
                            max_quality = max(phred_scores)
                            print(f"Quality stats: avg={avg_quality:.1f}, min={min_quality}, max={max_quality}")
                    
                    # Show masking statistics if threshold is used
                    if args.quality_threshold is not None:
                        phred_scores = [ord(q) - 33 for q in quality_string]
                        low_quality_count = sum(1 for q in phred_scores if q < args.quality_threshold)
                        masking_percentage = (low_quality_count / len(sequence)) * 100
                        print(f"Masked:   {low_quality_count}/{len(sequence)} bases ({masking_percentage:.1f}%)")
                    
                    # Show average quality if filtering is enabled
                    if args.avg_quality_threshold is not None:
                        avg_quality = calculate_average_quality(quality_string)
                        print(f"Avg Quality: {avg_quality:.1f} (threshold: {args.avg_quality_threshold})")
                
                # QUIET MODE: Only show sequence (and quality if requested), no formatting
                else:
                    print(display_sequence)
                    if args.show_quality:
                        if args.quality_threshold is not None:
                            formatted_quality = format_quality_scores(quality_string, args.quality_threshold)
                            print(formatted_quality)
                        else:
                            print(quality_string)
                
                # If looking for specific read ID, stop after finding it
                if args.read_id and args.read_id in record.id:
                    break
    
    except Exception as e:
        if args.verbose:
            logger.error(f"Error reading FASTQ file: {e}")
        sys.exit(1)
    
    finally:
        # Clean up temporary file if it was created
        if temp_file and os.path.exists(temp_file):
            os.unlink(temp_file)
    
    if args.verbose:
        print(f"\n{'='*80}")
        if args.random:
            print(f"Displayed {displayed_count} reads (target: {args.nreads})")
        else:
            print(f"Displayed {displayed_count} reads from {original_fastq}")
        
        if args.avg_quality_threshold is not None and not args.random:
            print(f"Filtered out {filtered_count} reads with average quality < {args.avg_quality_threshold}")
        
        if args.quality_threshold is not None:
            print(f"Used quality threshold: {args.quality_threshold} (bases below this are masked)")

def adaptive_sample_with_quality_filter(fastq_file, args):
    """Adaptively sample reads until we get enough high-quality ones"""
    logger = logging.getLogger(__name__)
    import tempfile
    import subprocess
    
    target_reads = args.nreads
    max_attempts = 5
    attempt = 0
    sample_multiplier = 3  # Start with 3x
    total_sampled = 0
    
    # Create final output file
    final_temp = tempfile.NamedTemporaryFile(mode='w+', suffix='.fastq', delete=False)
    final_path = final_temp.name
    final_temp.close()
    
    try:
        while attempt < max_attempts:
            attempt += 1
            
            # Calculate sample size for this attempt
            if attempt == 1:
                sample_size = target_reads * sample_multiplier
            else:
                # Increase sample size based on previous success rate
                sample_size = target_reads * sample_multiplier * attempt
            
            if args.verbose:
                logger.info(f"Attempt {attempt}: Sampling {sample_size} reads to find {target_reads} high-quality reads")
            
            # Sample reads
            temp_sample = sample_fastq_with_seqtk(fastq_file, sample_size, args.seed + attempt - 1)
            
            # Filter by quality and add to final file
            good_reads = 0
            with open_fastq(temp_sample) as sample_handle, open(final_path, 'a') as final_handle:
                for record in SeqIO.parse(sample_handle, "fastq"):
                    # Get quality string
                    quality_scores = record.letter_annotations.get("phred_quality", None)
                    if quality_scores is None:
                        quality_string = record.format("fastq").split('\n')[3]
                    else:
                        quality_string = ''.join(chr(q + 33) for q in quality_scores)
                    
                    # Check quality
                    avg_quality = calculate_average_quality(quality_string)
                    if avg_quality >= args.avg_quality_threshold:
                        # Write this read to final file
                        SeqIO.write(record, final_handle, "fastq")
                        good_reads += 1
                        
                        # Stop if we have enough
                        if good_reads >= target_reads:
                            break
            
            # Clean up temp sample file
            if os.path.exists(temp_sample):
                os.unlink(temp_sample)
            
            total_sampled += sample_size
            
            if args.verbose:
                logger.info(f"Found {good_reads} high-quality reads in attempt {attempt}")
            
            # Check if we have enough reads
            if good_reads >= target_reads:
                if args.verbose:
                    logger.info(f"Successfully found {good_reads} high-quality reads after {attempt} attempts (sampled {total_sampled} total)")
                break
            
            # If this is the last attempt, warn but continue
            if attempt == max_attempts:
                if args.verbose:
                    logger.warning(f"Could only find {good_reads} high-quality reads after {max_attempts} attempts")
                break
        
        return final_path
        
    except Exception as e:
        # Clean up on error
        if os.path.exists(final_path):
            os.unlink(final_path)
        raise e
       
def main():
    logger = setup_logging()
    
    try:
        fastq_file, args = get_parameters(logger)
        
        if args.verbose:
            logger.info(f"Input FASTQ: {fastq_file}")
        
        # View FASTQ sequences
        view_fastq_sequences(fastq_file, args)
        
        if args.verbose:
            logger.info("FASTQ viewing completed successfully!")
        
    except Exception as e:
        if args.verbose:
            logger.error(f"Error during FASTQ viewing: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
