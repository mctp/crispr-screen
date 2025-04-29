import os
import sys
import argparse
import subprocess
import logging
# project-specific imports
import utilities.parse_config as pc

def mageck_test(config_file='', count_matrix_file='', control_gene_id='', control_sample_min_count='', gene_sgRNA_min='', comparisons=None, normalization_method='median', remove_zero_method='control', output_dir=''):

    if config_file:
        config = pc.parse_config(config_file)
        output_dir = config.get('output_dir')
        experiment_name = config.get('experiment_name')
        if config.get('count_matrix_file'):
            count_matrix_file = config.get('count_matrix_file')
        else:
            if experiment_name is not None and experiment_name != '':
                count_matrix_file = os.path.join(output_dir, f"{experiment_name}_sgrna_count_matrix.txt")
            else:
                count_matrix_file = os.path.join(output_dir, "sgrna_count_matrix.txt")
        control_gene_id = config.get('control_gene_id')
        core_essential_genes_file = config.get('core_essential_genes_file')
        non_essential_genes_file = config.get('non_essential_genes_file')
        control_sample = config.get('control_sample')
        comparisons = config.get('comparisons')
    else:
        logging.info("No config file provided. Using command line arguments.")

    # defaults
    if config:
        normalization_method = normalization_method if normalization_method else config.get('normalization_method', 'median')
        remove_zero_method = remove_zero_method if remove_zero_method else config.get('remove_zero_method', 'control')
    else:
        normalization_method = normalization_method if normalization_method else 'median'
        remove_zero_method = remove_zero_method if remove_zero_method else 'control'

    logging.info(f"Normalization method: {normalization_method}")

    if not os.path.isfile(count_matrix_file) or os.path.getsize(count_matrix_file) == 0:
        logging.error(f"Error: Count matrix file '{count_matrix_file}' does not exist or is empty.")
        sys.exit(1)
    
    logging.info(f"Experiment: {experiment_name}")
    for comparison in comparisons:
        logging.info(f"   Comparison: {comparison}")
        a, b = comparison.split(':')
        
        with open(count_matrix_file, 'r') as f:
            header = f.readline().strip()
        
        if a not in header or b not in header:
            logging.error(f"Error: One or both samples ({a}, {b}) not found in the counts matrix header. Skipping comparison.")
            continue

        if normalization_method not in ['median', 'total', 'control']:
            logging.error(f"Error: Invalid normalization method '{normalization_method}'. Valid options are: median, total.")
            continue

        if remove_zero_method not in ['none', 'control', 'treatment', 'both', 'any']:
            logging.error(f"Error: Invalid remove_zero_method '{remove_zero_method}'. Valid options are: none, control, treatment, both, any.")
            continue

        comparison_name = f"{a}_vs_{b}"
        comparison_dir = os.path.join(output_dir, comparison_name)
        counter = 1
        while os.path.isdir(comparison_dir):
            comparison_dir = os.path.join(output_dir, f"{comparison_name}_{counter:02d}")
            counter += 1
        comparison_prefix = os.path.join(comparison_dir, comparison_name)
        os.makedirs(comparison_dir, exist_ok=True)

        # make control gene file if it doesn't exist
        control_gene_file = "control_gene.txt"
        if not control_gene_file or not os.path.isfile(control_gene_file) or os.path.getsize(control_gene_file) == 0:
            with open(control_gene_file, 'w') as f:
                f.write(f"{control_gene_id}\n")

        subprocess.run([
            'mageck', 'test',
            '-k', count_matrix_file,
            '-t', a,
            '-c', b,
            '--control-gene', control_gene_file,
            '-n', comparison_prefix,
            '--norm-method', normalization_method,
            '--remove-zero', remove_zero_method,
            '--normcounts-to-file'
        ])

        summary_file = os.path.join(comparison_dir, 'pipeline_variable_summary.txt')
        with open(summary_file, 'w') as f:
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

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Mageck Analysis')
    parser.add_argument('--config', help='Path to the config file')
    parser.add_argument('--mode', required=True, help='Mageck mode')
    parser.add_argument('--count_matrix_file', default='', help='Count matrix file')
    parser.add_argument('--experiment_name', required=True, help='Experiment name')
    parser.add_argument('--core_essential_genes_file', required=True, help='Core essential genes file')
    parser.add_argument('--non_essential_genes_file', required=True, help='Non essential genes file')
    parser.add_argument('--control_gene_id', default='', help='Control gene ID')
    parser.add_argument('--control_sample', required=True, help='Control sample')
    parser.add_argument('--control_sample_min_count', required=True, help='Control sample min count')
    parser.add_argument('--gene_sgRNA_min', required=True, help='Min sgRNAs per gene')
    parser.add_argument('--comparisons', nargs='+', help='Comparisons to process')
    parser.add_argument('--normalization_method', default='median', help='Normalization method')
    parser.add_argument('--remove_zero_method', default='control', help='Remove zero method')
    parser.add_argument('--output_dir', required=True, help='Output directory')
    args = parser.parse_args()

    mageck_test(
        config_file=args.config,
        count_matrix_file=args.count_matrix_file,
        control_gene_id=args.control_gene_id,
        control_sample_min_count=args.control_sample_min_count,
        gene_sgRNA_min=args.gene_sgRNA_min,
        comparisons=args.comparisons,
        normalization_method=args.normalization_method,
        remove_zero_method=args.remove_zero_method,
        output_dir=args.output_dir
        )
