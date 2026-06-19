import os
import sys
import argparse
import subprocess
import logging
import re
import pandas as pd
import numpy as np

# project-specific imports
import utilities.parse_config as pc

def sanitize_filename(name):
    """Sanitize analysis name for use as filename"""
    # Replace spaces with underscores and remove non-alphanumeric characters except underscore and hyphen
    sanitized = re.sub(r'[^\w\-]', '_', name)
    # Remove multiple consecutive underscores
    sanitized = re.sub(r'_+', '_', sanitized)
    # Remove leading/trailing underscores
    sanitized = sanitized.strip('_')
    return sanitized

def create_design_matrix(design_matrix_config, analysis_dir, experiment_name, analysis_name):
    """Create design matrix file for MAGeCK MLE from config"""
    
    # Create design matrix filename with analysis name
    safe_analysis_name = sanitize_filename(analysis_name)
    design_matrix_file = os.path.join(analysis_dir, f"{experiment_name}_{safe_analysis_name}.design_matrix.txt")
    
    with open(design_matrix_file, 'w') as f:
        samples = design_matrix_config['samples']
        conditions = design_matrix_config['conditions']
        
        # Header line: "Samples" + condition names
        condition_names = [condition['name'] for condition in conditions]
        f.write("Samples\t" + "\t".join(condition_names) + "\n")
        
        # Data rows: one row per sample
        for i, sample in enumerate(samples):
            # Start with sample name
            row = [sample]
            
            # Add value for each condition
            for condition in conditions:
                condition_values = condition['values']
                if len(condition_values) != len(samples):
                    raise ValueError(f"Design matrix condition '{condition['name']}' has {len(condition_values)} values but {len(samples)} samples")
                row.append(str(condition_values[i]))
            
            f.write("\t".join(row) + "\n")
    
    logging.info(f"Created design matrix file: {design_matrix_file}")
    return design_matrix_file

def create_control_gene_file(control_gene_id, analysis_dir):
    """Create control gene file for MAGeCK RRA and MLE"""
    control_gene_file = os.path.join(analysis_dir, "control_gene_list.txt")
    with open(control_gene_file, 'w') as f:
        f.write(f"{control_gene_id}\n")
    return control_gene_file

def validate_rra_parameters(normalization_method, remove_zero_method):
    """Validate RRA-specific parameters"""
    if normalization_method is not None and normalization_method not in ['median', 'total', 'control']:
        raise ValueError(f"Invalid normalization method '{normalization_method}'. Valid options are: median, total, control, or None (for MAGeCK default).")
    
    if remove_zero_method not in ['none', 'control', 'treatment', 'both', 'any']:
        raise ValueError(f"Invalid remove_zero_method '{remove_zero_method}'. Valid options are: none, control, treatment, both, any.")

def check_samples_in_matrix(samples, count_matrix_file):
    """Check if samples exist in count matrix header"""
    with open(count_matrix_file, 'r') as f:
        header = f.readline().strip()
    
    missing_samples = [sample for sample in samples if sample not in header]
    return missing_samples

def check_rra_output_exists(comparison_dir, comparison_name):
    """Check if RRA output files already exist"""
    gene_summary_file = os.path.join(comparison_dir, f"{comparison_name}.gene_summary.txt")
    sgrna_summary_file = os.path.join(comparison_dir, f"{comparison_name}.sgrna_summary.txt")
    
    return os.path.exists(gene_summary_file) and os.path.exists(sgrna_summary_file)

def check_mle_output_exists(mle_dir, mle_prefix):
    """Check if MLE output files already exist"""
    gene_summary_file = f"{mle_prefix}.gene_summary.txt"
    sgrna_summary_file = f"{mle_prefix}.sgrna_summary.txt"
    
    return os.path.exists(gene_summary_file) and os.path.exists(sgrna_summary_file)

def filter_sgrnas_for_rra(config, count_matrix_file, treatment_sample, control_sample):
    """Filter sgRNAs for RRA analysis based on rra_min_count, rra_min_cpm, and rra_min_count_method.
    If rra_min_cpm is set, use CPM matrix (compute if needed, or load if exists).
    If both are set, prefer count filtering for backward compatibility.
    """
    rra_min_count = config.get('rra_min_count')
    rra_min_cpm = config.get('rra_min_cpm')
    rra_min_count_method = config.get('rra_min_count_method', 'any')

    # If neither filtering parameter is specified, return None (no filtering)
    if not rra_min_count and not rra_min_cpm:
        return None

    # Prefer count filtering if both are set
    use_cpm = False
    threshold = None
    matrix_file = count_matrix_file

    if rra_min_count:
        threshold = rra_min_count
        use_cpm = False
        matrix_file = count_matrix_file
    elif rra_min_cpm:
        threshold = rra_min_cpm
        use_cpm = True
        # Determine CPM matrix file location
        cpm_matrix_file = config.get('cpm_matrix_file')
        if not cpm_matrix_file:
            output_dir = config.get('output_dir', 'output')
            experiment_name = config.get('experiment_name', '')
            if experiment_name:
                cpm_matrix_file = os.path.join(output_dir, f"{experiment_name}_sgrna_cpm_matrix.txt")
            else:
                cpm_matrix_file = os.path.join(output_dir, "sgrna_cpm_matrix.txt")
        matrix_file = cpm_matrix_file

        # Compute CPM matrix if it doesn't exist
        if not os.path.exists(matrix_file):
            import cpm_matrix as cpm
            logging.info(f"CPM matrix file {matrix_file} not found. Computing CPM matrix...")
            cpm.count_and_cpm_matrix(config_file=config.get('config_file', 'config.yaml'))
            if not os.path.exists(matrix_file):
                logging.error(f"Failed to generate CPM matrix file: {matrix_file}")
                return None

    # Load the appropriate matrix
    try:
        df = pd.read_csv(matrix_file, sep='\t', index_col=0)
    except Exception as e:
        logging.error(f"Error loading matrix file {matrix_file}: {e}")
        return None

    # Check that specified samples exist in the matrix
    missing_samples = [s for s in [treatment_sample, control_sample] if s not in df.columns]
    if missing_samples:
        logging.error(f"Specified samples not found in matrix: {missing_samples}")
        return None

    # Apply filtering based on method
    if rra_min_count_method == 'any':
        mask = (df[treatment_sample] >= threshold) | (df[control_sample] >= threshold)
    elif rra_min_count_method == 'all':
        mask = (df[treatment_sample] >= threshold) & (df[control_sample] >= threshold)
    elif rra_min_count_method == 'control':
        mask = df[control_sample] >= threshold
    else:
        logging.error(f"Invalid rra_min_count_method '{rra_min_count_method}'. Valid options are: any, all, control")
        return None

    passing_sgrnas = df[mask].index.tolist()
    filter_type = "CPM" if use_cpm else "count"
    logging.info(f"RRA sgRNA filtering: {len(passing_sgrnas)}/{len(df)} sgRNAs passed "
                 f"{filter_type} threshold of {threshold} using method '{rra_min_count_method}' "
                 f"for samples {treatment_sample}, {control_sample}")

    return passing_sgrnas

def run_mageck_rra(config, comparisons, count_matrix_file, output_dir, experiment_name, force_rerun=False):
    """Run MAGeCK RRA analysis for all comparisons"""
    
    logging.info("Running MAGeCK RRA analysis")
    logging.info(f"Total RRA comparisons to process: {len(comparisons)}")
    
    # Get RRA-specific configuration with defaults
    control_gene_id = config.get('control_gene_id', '')
    core_essential_genes_file = config.get('core_essential_genes_file', '')
    non_essential_genes_file = config.get('non_essential_genes_file', '')
    control_sample = config.get('control_sample', '')
    control_sample_min_count = config.get('control_sample_min_count', '')
    gene_sgRNA_min = config.get('gene_sgRNA_min', '')
    normalization_method = config.get('normalization_method', 'median')
    remove_zero_method = config.get('remove_zero_method', 'control')
    rra_min_good_sgrna = config.get('rra_min_good_sgrna')
    
    # Handle None/empty normalization method
    if normalization_method in [None, '', 'null', 'none']:
        normalization_method = None
        logging.info("Normalization method: None (will use MAGeCK default)")
    else:
        logging.info(f"Normalization method: {normalization_method}")

    # Validate parameters once
    try:
        validate_rra_parameters(normalization_method, remove_zero_method)
    except ValueError as e:
        logging.error(f"Error: {e}")
        return
    
    processed_count = 0
    skipped_count = 0
    failed_count = 0
    
    for comparison_item in comparisons:
        processed_count += 1
        
        # Handle both legacy and new comparison formats
        if isinstance(comparison_item, str):
            # Legacy format: "treatment:control"
            comparison = comparison_item
            suffix = ""
        elif isinstance(comparison_item, dict):
            # New format: {"comparison": "treatment:control", "suffix": "_30cpm"}
            comparison = comparison_item.get('comparison', '')
            suffix = comparison_item.get('suffix', '')
            if not comparison:
                logging.error(f"Invalid comparison format: missing 'comparison' key in {comparison_item}. Skipping.")
                failed_count += 1
                continue
        else:
            logging.error(f"Invalid comparison format: {comparison_item}. Expected string or dict. Skipping.")
            failed_count += 1
            continue
            
        logging.info(f"[{processed_count}/{len(comparisons)}] Processing comparison: {comparison}")
        
        try:
            a, b = comparison.split(':')
        except ValueError:
            logging.error(f"Invalid comparison format '{comparison}'. Expected 'sample1:sample2'. Skipping.")
            failed_count += 1
            continue
        
        # Check if samples exist in matrix
        missing_samples = check_samples_in_matrix([a, b], count_matrix_file)
        if missing_samples:
            logging.error(f"Samples not found in count matrix: {missing_samples}. Skipping comparison {comparison}.")
            failed_count += 1
            continue

        # Create comparison directory name with new naming scheme
        comparison_name = f"{a}_vs_{b}"
        dir_name = f"{experiment_name}_rra_{comparison_name}{suffix}"
        comparison_dir = os.path.join(output_dir, dir_name)
        
        # Handle directory conflicts by appending counter
        counter = 1
        original_dir = comparison_dir
        while os.path.isdir(comparison_dir) and not check_rra_output_exists(comparison_dir, comparison_name):
            comparison_dir = f"{original_dir}_{counter:02d}"
            counter += 1
        
        # Check if analysis already exists and should be skipped
        if not force_rerun and os.path.isdir(comparison_dir) and check_rra_output_exists(comparison_dir, comparison_name):
            logging.info(f"[{processed_count}/{len(comparisons)}] RRA analysis for {comparison} already exists in {comparison_dir}. Skipping. Use --force to rerun.")
            skipped_count += 1
            continue
        
        comparison_prefix = os.path.join(comparison_dir, comparison_name)
        os.makedirs(comparison_dir, exist_ok=True)

        # Filter sgRNAs if filtering parameters are specified
        filtered_sgrnas = filter_sgrnas_for_rra(config, count_matrix_file, a, b)
        
        # Create filtered count matrix if filtering was applied
        analysis_count_matrix_file = count_matrix_file
        if filtered_sgrnas is not None:
            # Load original count matrix
            count_df = pd.read_csv(count_matrix_file, sep='\t', index_col=0)
            
            # Filter to only include passing sgRNAs
            filtered_count_df = count_df.loc[filtered_sgrnas]
            
            # Save filtered matrix
            filtered_count_file = os.path.join(comparison_dir, f"{comparison_name}_filtered_count_matrix.txt")
            filtered_count_df.to_csv(filtered_count_file, sep='\t')
            
            # Use filtered matrix for analysis
            analysis_count_matrix_file = filtered_count_file
            logging.info(f"Using filtered count matrix: {filtered_count_file}")

        # Create control gene file (only if control_gene_id is set)
        control_gene_file = None
        if control_gene_id:
            control_gene_file = create_control_gene_file(control_gene_id, comparison_dir)

        # Run MAGeCK test
        cmd = [
            'mageck', 'test',
            '-k', analysis_count_matrix_file,
            '-t', a,
            '-c', b,
            '-n', comparison_prefix,
            '--remove-zero', remove_zero_method,
            '--normcounts-to-file',
            '--keep-tmp'
        ]
        # Add control gene file only if present
        if control_gene_file:
            cmd.extend(['--control-gene', control_gene_file])
        # Add normalization method only if specified
        if normalization_method is not None:
            cmd.extend(['--norm-method', normalization_method])
        # Add additional-rra-parameters if rra_min_good_sgrna is set and not empty
        if rra_min_good_sgrna is not None and str(rra_min_good_sgrna).strip() != '':
            cmd.extend(['--additional-rra-parameters', f'--min-number-goodsgrna {rra_min_good_sgrna}'])

        logging.info(f"[{processed_count}/{len(comparisons)}] Running MAGeCK test for {comparison}")
        print(f"Running MAGeCK RRA analysis {processed_count}/{len(comparisons)}: {comparison}", flush=True)
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error(f"[{processed_count}/{len(comparisons)}] MAGeCK test failed for comparison {comparison} with return code {result.returncode}")
            logging.error(f"STDERR: {result.stderr}")
            failed_count += 1
            continue

        logging.info(f"[{processed_count}/{len(comparisons)}] MAGeCK RRA completed successfully for comparison {comparison}")
        print(f"✓ Completed MAGeCK RRA analysis {processed_count}/{len(comparisons)}: {comparison}", flush=True)

        # Write summary file
        write_rra_summary(comparison_dir, comparison_name, experiment_name, analysis_count_matrix_file,
                         core_essential_genes_file, non_essential_genes_file, control_gene_id,
                         control_sample, control_sample_min_count, gene_sgRNA_min, comparison,
                         normalization_method, remove_zero_method)
    
    # Summary
    successful_count = processed_count - skipped_count - failed_count
    logging.info(f"RRA Analysis Summary: {successful_count} completed, {skipped_count} skipped, {failed_count} failed")
    print(f"RRA Analysis Summary: {successful_count} completed, {skipped_count} skipped, {failed_count} failed", flush=True)

def write_rra_summary(comparison_dir, comparison_name, experiment_name, count_matrix_file,
                     core_essential_genes_file, non_essential_genes_file, control_gene_id,
                     control_sample, control_sample_min_count, gene_sgRNA_min, comparison,
                     normalization_method, remove_zero_method):
    """Write RRA analysis summary file"""
    summary_file = os.path.join(comparison_dir, 'pipeline_variable_summary.txt')
    with open(summary_file, 'w') as f:
        f.write(f"method: rra\n")
        f.write(f"experiment name: {experiment_name}\n")
        f.write(f"count matrix file: {count_matrix_file}\n")
        f.write(f"core essential genes file: {core_essential_genes_file}\n")
        f.write(f"non essential genes file: {non_essential_genes_file}\n")
        f.write(f"control gene ID: {control_gene_id}\n")
        f.write(f"control sample: {control_sample}\n")
        f.write(f"control sample min count: {control_sample_min_count}\n")
        f.write(f"min sgRNAs per gene: {gene_sgRNA_min}\n")
        f.write(f"comparison: {comparison}\n")
        f.write(f"normalization method: {normalization_method}\n")
        f.write(f"remove zero method: {remove_zero_method}\n")

def run_mageck_mle(config, design_matrices_config, count_matrix_file, cpm_matrix_file, output_dir, experiment_name, force_rerun=False):
    """Run MAGeCK MLE analysis for all design matrices"""
    
    logging.info("Running MAGeCK MLE analysis")
    logging.info(f"Total MLE analyses to process: {len(design_matrices_config)}")
    
    # Get MLE-specific configuration with defaults
    core_essential_genes_file = config.get('core_essential_genes_file', '')
    non_essential_genes_file = config.get('non_essential_genes_file', '')
    
    # Handle MLE threads with fallback to ncpu, then default
    mle_threads = config.get('mle_threads')
    if not mle_threads:  # If empty, None, or 0
        mle_threads = config.get('ncpu', 32)  # Fall back to ncpu, then default 32
    
    mle_permutations = config.get('mle_permutations', 10)
    
    # Use control_gene_id from config for MLE control gene
    mle_control_gene = config.get('control_gene_id', '')
    
    processed_count = 0
    skipped_count = 0
    failed_count = 0
    
    # Process each design matrix
    for i, design_matrix_config in enumerate(design_matrices_config):
        processed_count += 1
        analysis_name = design_matrix_config.get('analysis_name', f'analysis_{i+1}')
        safe_analysis_name = sanitize_filename(analysis_name)
        
        logging.info(f"[{processed_count}/{len(design_matrices_config)}] Processing analysis: {analysis_name}")
        
        try:
            # Create MLE output directory for this analysis
            mle_dir = os.path.join(output_dir, f"{experiment_name}_mle_{safe_analysis_name}")
            mle_prefix = os.path.join(mle_dir, f"{experiment_name}_{safe_analysis_name}")
            
            # Check if analysis already exists and should be skipped
            if not force_rerun and os.path.isdir(mle_dir) and check_mle_output_exists(mle_dir, mle_prefix):
                logging.info(f"[{processed_count}/{len(design_matrices_config)}] MLE analysis '{analysis_name}' already exists in {mle_dir}. Skipping. Use --force to rerun.")
                print(f"Skipping MAGeCK MLE analysis {processed_count}/{len(design_matrices_config)}: {analysis_name} (already exists)", flush=True)
                skipped_count += 1
                continue
            
            os.makedirs(mle_dir, exist_ok=True)
            
            # Filter sgRNAs if filtering parameters are specified
            filtered_sgrnas = filter_sgrnas_by_threshold(design_matrix_config, count_matrix_file, cpm_matrix_file)
            
            # Create filtered count matrix if filtering was applied
            analysis_count_matrix_file = count_matrix_file
            if filtered_sgrnas is not None:
                # Load original count matrix
                count_df = pd.read_csv(count_matrix_file, sep='\t', index_col=0)
                
                # Filter to only include passing sgRNAs
                filtered_count_df = count_df.loc[filtered_sgrnas]
                
                # Save filtered matrix
                filtered_count_file = os.path.join(mle_dir, f"{experiment_name}_{safe_analysis_name}_filtered_count_matrix.txt")
                filtered_count_df.to_csv(filtered_count_file, sep='\t')
                
                # Use filtered matrix for analysis
                analysis_count_matrix_file = filtered_count_file
                logging.info(f"Using filtered count matrix: {filtered_count_file}")
            
            # Create design matrix from config (now in analysis directory)
            design_matrix_file = create_design_matrix(design_matrix_config, mle_dir, experiment_name, analysis_name)
            
            # Create control gene file
            control_gene_file = create_control_gene_file(mle_control_gene, mle_dir)
            
            # Build and run MAGeCK MLE command
            cmd = [
                'mageck', 'mle',
                '-k', analysis_count_matrix_file,
                '-d', design_matrix_file,
                '-n', mle_prefix,
                '--control-gene', control_gene_file,
                '--threads', str(mle_threads),
                '--update-efficiency',
                '--permutation-round', str(mle_permutations)
            ]
            
            logging.info(f"[{processed_count}/{len(design_matrices_config)}] Running MAGeCK MLE for analysis '{analysis_name}'")
            print(f"Running MAGeCK MLE analysis {processed_count}/{len(design_matrices_config)}: {analysis_name}", flush=True)
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                logging.error(f"[{processed_count}/{len(design_matrices_config)}] MAGeCK MLE failed for analysis '{analysis_name}' with return code {result.returncode}")
                logging.error(f"STDERR: {result.stderr}")
                failed_count += 1
                continue  # Continue with next analysis instead of exiting
            else:
                logging.info(f"[{processed_count}/{len(design_matrices_config)}] MAGeCK MLE completed successfully for analysis '{analysis_name}'")
                logging.info(f"Output files written to: {mle_dir}")
                print(f"✓ Completed MAGeCK MLE analysis {processed_count}/{len(design_matrices_config)}: {analysis_name}", flush=True)
            
            # Write summary file
            write_mle_summary(mle_dir, analysis_name, experiment_name, analysis_count_matrix_file,
                             design_matrix_file, core_essential_genes_file, non_essential_genes_file,
                             mle_control_gene, mle_threads, mle_permutations, design_matrix_config)
                             
        except Exception as e:
            logging.error(f"[{processed_count}/{len(design_matrices_config)}] Error processing MLE analysis '{analysis_name}': {e}")
            print(f"Failed MAGeCK MLE analysis {processed_count}/{len(design_matrices_config)}: {analysis_name} (error: {e})", flush=True)
            failed_count += 1
            continue
    
    # Summary
    successful_count = processed_count - skipped_count - failed_count
    logging.info(f"MLE Analysis Summary: {successful_count} completed, {skipped_count} skipped, {failed_count} failed")
    print(f"MLE Analysis Summary: {successful_count} completed, {skipped_count} skipped, {failed_count} failed", flush=True)

def write_mle_summary(mle_dir, analysis_name, experiment_name, count_matrix_file,
                     design_matrix_file, core_essential_genes_file, non_essential_genes_file,
                     mle_control_gene, mle_threads, mle_permutations, design_matrix_config):
    """Write MLE analysis summary file"""
    summary_file = os.path.join(mle_dir, 'pipeline_variable_summary.txt')
    with open(summary_file, 'w') as f:
        f.write(f"method: mle\n")
        f.write(f"analysis_name: {analysis_name}\n")
        f.write(f"experiment name: {experiment_name}\n")
        f.write(f"count matrix file: {count_matrix_file}\n")
        f.write(f"design matrix file: {design_matrix_file}\n")
        f.write(f"core essential genes file: {core_essential_genes_file}\n")
        f.write(f"non essential genes file: {non_essential_genes_file}\n")
        f.write(f"control gene ID: {mle_control_gene}\n")
        f.write(f"threads: {mle_threads}\n")
        f.write(f"permutation rounds: {mle_permutations}\n")
        f.write(f"samples: {', '.join(design_matrix_config['samples'])}\n")
        f.write(f"conditions: {', '.join([c['name'] for c in design_matrix_config['conditions']])}\n")

def get_design_matrices_config(config):
    """Get design matrices configuration with legacy fallback"""
    design_matrices_config = config.get('design_matrices')
    
    # Fallback to legacy single design matrix
    if not design_matrices_config:
        single_design_matrix = config.get('design_matrix')
        if single_design_matrix:
            logging.warning("Using legacy 'design_matrix' config. Consider migrating to 'design_matrices' for multiple analysis support.")
            design_matrices_config = [{
                'analysis_name': 'default',
                'samples': single_design_matrix['samples'],
                'conditions': single_design_matrix['conditions']
            }]
    
    return design_matrices_config

def filter_sgrnas_by_threshold(design_matrix, count_matrix_file, cpm_matrix_file=None):
    """Filter sgRNAs based on min_count, min_cpm, and min_sample parameters"""
    sgrna_min_count = design_matrix.get('sgrna_min_count')
    sgrna_min_cpm = design_matrix.get('sgrna_min_cpm')
    sgrna_min_sample = design_matrix.get('sgrna_min_sample')
    
    # If no filtering parameters specified, return None (no filtering)
    if not sgrna_min_count and not sgrna_min_cpm:
        return None
    
    # Determine which threshold to use (prefer min_count if both are specified)
    use_cpm = False
    threshold = None
    matrix_file = count_matrix_file
    
    if sgrna_min_count:
        threshold = sgrna_min_count
        use_cpm = False
        matrix_file = count_matrix_file
    elif sgrna_min_cpm:
        threshold = sgrna_min_cpm
        use_cpm = True
        matrix_file = cpm_matrix_file
        if not matrix_file or not os.path.exists(matrix_file):
            logging.error(f"CPM matrix file not found: {matrix_file}")
            return None
    
    # Load the appropriate matrix
    try:
        df = pd.read_csv(matrix_file, sep='\t', index_col=0)
    except Exception as e:
        logging.error(f"Error loading matrix file {matrix_file}: {e}")
        return None
    
    # Determine which samples to check
    samples_to_check = []
    if sgrna_min_sample:
        if isinstance(sgrna_min_sample, str):
            samples_to_check = [sgrna_min_sample]
        elif isinstance(sgrna_min_sample, list):
            samples_to_check = sgrna_min_sample
        else:
            samples_to_check = [str(sgrna_min_sample)]
    else:
        # If no specific sample specified, use all samples
        samples_to_check = df.columns.tolist()
        logging.warning("No sgrna_min_sample specified - filtering based on ALL samples. "
                       "This is unusual and could cause many sgRNAs to be rejected.")
    
    # Check that specified samples exist in the matrix
    missing_samples = [s for s in samples_to_check if s not in df.columns]
    if missing_samples:
        logging.error(f"Specified samples not found in matrix: {missing_samples}")
        return None
    
    # Apply filtering
    if len(samples_to_check) == 1:
        # Single sample - check if it meets threshold
        mask = df[samples_to_check[0]] >= threshold
    else:
        # Multiple samples - all must meet threshold
        mask = (df[samples_to_check] >= threshold).all(axis=1)
    
    passing_sgrnas = df[mask].index.tolist()
    
    filter_type = "CPM" if use_cpm else "count"
    logging.info(f"sgRNA filtering: {len(passing_sgrnas)}/{len(df)} sgRNAs passed "
                f"{filter_type} threshold of {threshold} in samples {samples_to_check}")
    
    return passing_sgrnas

def mageck_test(config_file='', count_matrix_file='', control_gene_id='', control_sample_min_count='', 
                gene_sgRNA_min='', comparisons=None, normalization_method='median', 
                remove_zero_method='control', output_dir='', method='auto', force_rerun=False, metadata_file=None):

    if config_file:
        config = pc.parse_config(config_file)
        output_dir = config.get('output_dir')
        experiment_name = config.get('experiment_name')

        if config.get('count_matrix_file'):
            count_matrix_file = config.get('count_matrix_file')
        else:
            if experiment_name:
                count_matrix_file = os.path.join(output_dir, f"{experiment_name}_sgrna_count_matrix.txt")
            else:
                count_matrix_file = os.path.join(output_dir, "sgrna_count_matrix.txt")

        if config.get('cpm_matrix_file'):
            cpm_matrix_file = config.get('cpm_matrix_file')
        else:
            if experiment_name:
                cpm_matrix_file = os.path.join(output_dir, f"{experiment_name}_sgrna_cpm_matrix.txt")
            else:
                cpm_matrix_file = os.path.join(output_dir, "sgrna_cpm_matrix.txt")

        comparisons = config.get('comparisons')
    else:
        logging.info("No config file provided. Using command line arguments.")
        config = {}  # Empty config for command line mode

    # Normalize method parameter
    method = method.lower()
    if method == 'test':
        method = 'rra'

    # Determine which analyses to run based on configuration
    run_mle = False
    run_rra = False
    design_matrices_config = None
    
    # Check for MLE design matrices
    if config_file:
        design_matrices_config = get_design_matrices_config(config)
        if design_matrices_config:
            run_mle = True
            logging.info(f"Found {len(design_matrices_config)} design matrices to process")
            print(f"✓ MLE analysis available: {len(design_matrices_config)} design matrices found", flush=True)
        else:
            print("✗ MLE analysis not available: no design matrices found", flush=True)
    
    # Check for RRA comparisons
    if comparisons:
        run_rra = True
        logging.info(f"Found {len(comparisons)} RRA comparisons to process")
        print(f"✓ RRA analysis available: {len(comparisons)} comparisons found", flush=True)
    else:
        print("✗ RRA analysis not available: no comparisons found", flush=True)
    
    print(f"Available analyses: MLE={run_mle}, RRA={run_rra}", flush=True)
    
    # Only use config analysis_method if explicitly specified AND command line method is auto
    if method == 'auto':
        config_method = config.get('analysis_method')
        if config_method and config_method.strip():  # Only use if not empty/whitespace
            config_method = config_method.lower().strip()
            print(f"Config specifies analysis_method: {config_method}", flush=True)
            # Check if the config method is compatible with available analyses
            if config_method == 'mle' and not run_mle:
                logging.warning(f"Config specifies analysis_method '{config_method}' but no design matrices found. Will auto-detect instead.")
                print(f"⚠ Config analysis_method '{config_method}' incompatible, falling back to auto-detect", flush=True)
            elif config_method == 'rra' and not run_rra:
                logging.warning(f"Config specifies analysis_method '{config_method}' but no comparisons found. Will auto-detect instead.")
                print(f"⚠ Config analysis_method '{config_method}' incompatible, falling back to auto-detect", flush=True)
            elif config_method == 'both' and not (run_mle and run_rra):
                logging.warning(f"Config specifies analysis_method '{config_method}' but missing analyses. Will auto-detect instead.")
                print(f"⚠ Config analysis_method '{config_method}' incompatible, falling back to auto-detect", flush=True)
            else:
                method = config_method
                logging.info(f"Using analysis_method from config: {method}")
                print(f"Using analysis_method from config: {method}", flush=True)
        else:
            print("No analysis_method specified in config (or empty), will auto-detect", flush=True)
    
    # Auto-detect method based on what's configured
    if method == 'auto':
        print("Auto-detecting method based on available analyses...", flush=True)
        if run_mle and run_rra:
            method = 'both'
            logging.info("Auto-detected method: BOTH (MLE + RRA)")
            print("Auto-detected method: BOTH (will run MLE + RRA)", flush=True)
        elif run_mle:
            method = 'mle'
            logging.info("Auto-detected method: MLE")
            print("Auto-detected method: MLE (only design matrices found)", flush=True)
        elif run_rra:
            method = 'rra'
            logging.info("Auto-detected method: RRA")
            print("Auto-detected method: RRA (only comparisons found)", flush=True)
        else:
            logging.error("Error: No analyses configured. Please specify either 'design_matrices' (for MLE) or 'comparisons' (for RRA) in your config file.")
            print("✗ Error: No analyses configured in config file", flush=True)
            sys.exit(1)
    else:
        # Validate explicit method choice
        if method not in ['rra', 'mle', 'both']:
            logging.error(f"Error: Invalid method '{method}'. Valid options are: rra, mle, both, auto.")
            sys.exit(1)
        
        # Validate that we have the required configuration for the specified method
        if method == 'mle' and not run_mle:
            logging.error("Error: design_matrices configuration required for MLE method")
            sys.exit(1)
        elif method == 'rra' and not run_rra:
            logging.error("Error: comparisons configuration required for RRA method")
            sys.exit(1)
        elif method == 'both' and not (run_mle and run_rra):
            logging.error("Error: Both design_matrices and comparisons required for 'both' method")
            sys.exit(1)
    
    logging.info(f"Experiment: {experiment_name}")
    logging.info(f"MAGeCK analysis mode: {method.upper()}")
    print(f"Starting MAGeCK analysis mode: {method.upper()}", flush=True)

    # Validate count matrix file exists
    if not os.path.isfile(count_matrix_file) or os.path.getsize(count_matrix_file) == 0:
        logging.error(f"Error: Count matrix file '{count_matrix_file}' does not exist or is empty.")
        sys.exit(1)
    
    # Auto-detect metadata file with priority order: 1) command line, 2) config, 3) auto-detect YAML, 4) auto-detect TXT
    if metadata_file is None:
        # Check config file for metadata_file setting
        metadata_file = config.get("metadata_file", config.get("METADATA_FILE", None))
        
    if metadata_file is None:
        # Auto-detect: prefer YAML, fallback to TXT
        if os.path.exists("sample_metadata.yaml"):
            metadata_file = "sample_metadata.yaml"
            logging.info("Auto-detected metadata file: sample_metadata.yaml")
        elif os.path.exists("sample_metadata.txt"):
            metadata_file = "sample_metadata.txt"
            logging.info("Auto-detected metadata file: sample_metadata.txt")
        else:
            logging.error("Neither sample_metadata.yaml nor sample_metadata.txt found.")
            logging.error("Please create a metadata file or specify one in the config.")
            sys.exit(1)

    logging.info(f"Using metadata file: {metadata_file}")
    
    # Run analyses based on the determined method
    if method in ['mle', 'both'] and run_mle:
        print(f"Starting MLE analyses for {len(design_matrices_config)} design matrices...", flush=True)
        run_mageck_mle(config, design_matrices_config, count_matrix_file, cpm_matrix_file, output_dir, experiment_name, force_rerun)

    if method in ['rra', 'both'] and run_rra:
        print(f"Starting RRA analyses for {len(comparisons)} comparisons...", flush=True)
        run_mageck_rra(config, comparisons, count_matrix_file, output_dir, experiment_name, force_rerun)

    print("All MAGeCK analyses completed!", flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='MAGeCK Analysis')
    parser.add_argument('--config', help='Path to the config file (default: config.yaml if present)')
    parser.add_argument('--metadata', help='Path to the metadata file (default: auto-detect sample_metadata.yaml or sample_metadata.txt)')
    parser.add_argument('--method', default='auto', help='MAGeCK method: auto, rra, mle, both (default: auto - detect from config)')
    parser.add_argument('--mode', help='Mageck mode (deprecated, use --method)')
    parser.add_argument('--force', action='store_true', help='Force rerun of analyses even if output files exist')
    parser.add_argument('--count_matrix_file', default='', help='Count matrix file')
    
    # Make these arguments optional when config is provided
    parser.add_argument('--experiment_name', help='Experiment name (required if no config)')
    parser.add_argument('--core_essential_genes_file', help='Core essential genes file (required if no config)')
    parser.add_argument('--non_essential_genes_file', help='Non essential genes file (required if no config)')
    parser.add_argument('--control_gene_id', default='', help='Control gene ID')
    parser.add_argument('--control_sample', help='Control sample (required if no config)')
    parser.add_argument('--control_sample_min_count', help='Control sample min count (required if no config)')
    parser.add_argument('--gene_sgRNA_min', help='Min sgRNAs per gene (required if no config)')
    parser.add_argument('--comparisons', nargs='+', help='Comparisons to process (for RRA method)')
    parser.add_argument('--normalization_method', default='median', help='Normalization method (for RRA method)')
    parser.add_argument('--remove_zero_method', default='control', help='Remove zero method (for RRA method)')
    parser.add_argument('--output_dir', help='Output directory (required if no config)')
    args = parser.parse_args()

    # Determine config file to use
    config_file = args.config
    if not config_file:
        # Check for default config.yaml in working directory
        default_config = 'config.yaml'
        if os.path.exists(default_config):
            config_file = default_config
            logging.info(f"Using default config file: {default_config}")
    
    # Check if config file is provided/found and exists
    if config_file:
        if not os.path.exists(config_file):
            logging.error(f"Config file not found: {config_file}")
            sys.exit(1)
        
        # When config is provided, use it for everything
        logging.info(f"Using config file: {config_file}")
        
        # Handle deprecated --mode argument
        method = args.method
        if args.mode and args.method == 'auto':  # Only override if still default
            logging.warning("--mode is deprecated, use --method instead")
            method = args.mode
        
        mageck_test(
            config_file=config_file,
            method=method,
            force_rerun=args.force,
            metadata_file=args.metadata
        )
    else:
        # When no config is provided, require all parameters
        required_args = ['experiment_name', 'core_essential_genes_file', 'non_essential_genes_file', 
                        'control_sample', 'control_sample_min_count', 'gene_sgRNA_min', 'output_dir']
        
        missing_args = [arg for arg in required_args if getattr(args, arg) is None]
        if missing_args:
            logging.error(f"The following arguments are required when no config file is provided: {', '.join(['--' + arg for arg in missing_args])}")
            logging.error("Alternatively, place a 'config.yaml' file in the working directory or specify --config")
            parser.print_help()
            sys.exit(1)
        
        # Handle deprecated --mode argument
        method = args.method
        if args.mode and args.method == 'auto':  # Only override if still default
            logging.warning("--mode is deprecated, use --method instead")
            method = args.mode

        mageck_test(
            config_file='',
            count_matrix_file=args.count_matrix_file,
            control_gene_id=args.control_gene_id,
            control_sample_min_count=args.control_sample_min_count,
            gene_sgRNA_min=args.gene_sgRNA_min,
            comparisons=args.comparisons,
            normalization_method=args.normalization_method,
            remove_zero_method=args.remove_zero_method,
            output_dir=args.output_dir,
            method=method,
            force_rerun=args.force,
            metadata_file=args.metadata
        )

