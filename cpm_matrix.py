import argparse
import pandas as pd
import os
import logging
import sys
# project-specific imports
import utilities.utilities as util
import utilities.process_metadata as pm

def save_matrix_with_excel(df, output_prefix, matrix_type):
    """Save matrix to txt, csv, and xlsx formats with formatted Excel output"""
    # Save txt and csv as before
    df.to_csv(f"{output_prefix}sgrna_{matrix_type}_matrix.txt", sep='\t')
    df.to_csv(f"{output_prefix}sgrna_{matrix_type}_matrix.csv")
    
    # Save Excel with formatting
    excel_file = f"{output_prefix}sgrna_{matrix_type}_matrix.xlsx"
    with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
        df.to_excel(writer, sheet_name=f'{matrix_type.upper()}_Matrix', index=True)
        
        # Get the workbook and worksheet
        workbook = writer.book
        worksheet = writer.sheets[f'{matrix_type.upper()}_Matrix']
        
        # Create a bold font style
        from openpyxl.styles import Font, Border
        bold_font = Font(bold=True)
        normal_font = Font(bold=False)
        no_border = Border()
        
        # First, remove all borders and reset all fonts to normal
        for row in worksheet.iter_rows():
            for cell in row:
                cell.border = no_border
                cell.font = normal_font
        
        # Then apply bold formatting ONLY to the header row (row 1)
        for cell in worksheet[1]:
            cell.font = bold_font
        
        # Auto-size columns
        for column in worksheet.columns:
            max_length = 0
            column_letter = column[0].column_letter
            for cell in column:
                try:
                    if len(str(cell.value)) > max_length:
                        max_length = len(str(cell.value))
                except:
                    pass
            adjusted_width = min(max_length + 2, 50)  # Cap at 50 characters
            worksheet.column_dimensions[column_letter].width = adjusted_width
    
    logging.info(f"{matrix_type.capitalize()} matrix saved to {output_prefix}sgrna_{matrix_type}_matrix.txt, .csv, and .xlsx")

def count_matrix(count_file_paths):
    combined_df = pd.DataFrame()
    for file_path in count_file_paths:
        if os.path.exists(file_path):
            df = pd.read_csv(file_path, sep='\t', header=0)
            df.rename(columns={'sgRNA': 'sgrna_id', 'Gene': 'sgrna_target'}, inplace=True)
            df.set_index(['sgrna_id', 'sgrna_target'], inplace=True)
            if combined_df.empty:
                combined_df = df
            else:
                combined_df = combined_df.join(df, how='outer')
        else:
            logging.warning(f"File not found: {file_path}")
    combined_df.fillna(0, inplace=True)
    return combined_df

def cpm_matrix(count_matrix):
    cpm_matrix = count_matrix.div(count_matrix.sum(axis=0), axis=1) * 1e6
    return cpm_matrix


def count_and_cpm_matrix(config_file=None, metadata_file=None, input_fastq_dir=None, output_prefix=None, mode=None):
    
    logging.info(f"config_file: {config_file}")
    logging.info(f"metadata_file: {metadata_file}")
    logging.info(f"input_fastq_dir: {input_fastq_dir}")
    logging.info(f"output_prefix: {output_prefix}")
    logging.info(f"mode: {mode}")

    # Handle config file logic
    use_config = True
    if config_file is None:
        # Check for default config.yaml
        default_config = 'config.yaml'
        if os.path.exists(default_config):
            config_file = default_config
            logging.info(f"Using default config file: {config_file}")
        else:
            use_config = False
            logging.info("No config file specified and config.yaml not found. Running without config.")
    elif config_file.lower() == 'none':
        use_config = False
        config_file = None
        logging.info("Config explicitly disabled with --config none")

    if not use_config:
        # attempt to load without a config file
        if not metadata_file or not os.path.isfile(metadata_file):
            raise FileNotFoundError(f"Metadata file {metadata_file} does not exist or was not specified.")
        if not input_fastq_dir or not os.path.isdir(input_fastq_dir):
            raise NotADirectoryError(f"Input FASTQ directory {input_fastq_dir} does not exist, is not a directory, or was not specified.")
        if mode is None:
            raise ValueError("Mode is required if config is not provided.")
    else:
        logging.info(f"Loading config file {config_file}...")
        config = util.load_config(quiet=True, config_file=config_file)

    if use_config:
        # note: function params override config file
        metadata_file = metadata_file if metadata_file is not None else config.get('METADATA_FILE', config.get('metadata_file', None))
        input_fastq_dir = input_fastq_dir if input_fastq_dir is not None else config.get('FASTQ_DIR', config.get('fastq_dir', None))
        mode = mode if mode is not None else config.get('MODE', config.get('mode', None))
        output_dir = config.get('OUTPUT_DIR', config.get('output_dir', '.'))
        experiment_name = config.get('EXPERIMENT_NAME', config.get('experiment_name', 'experiment'))
    else:
        output_dir = '.'
        experiment_name = 'experiment'

    # Auto-detect metadata file if not specified by command line or config
    if metadata_file is None:
        if os.path.exists("sample_metadata.yaml"):
            metadata_file = "sample_metadata.yaml"
            logging.info("Auto-detected metadata file: sample_metadata.yaml")
        elif os.path.exists("sample_metadata.txt"):
            metadata_file = "sample_metadata.txt"
            logging.info("Auto-detected metadata file: sample_metadata.txt")
        else:
            raise FileNotFoundError("Neither sample_metadata.yaml nor sample_metadata.txt found. Please create a metadata file or specify one with --metadata.")

    output_prefix = output_prefix if output_prefix is not None else os.path.join(output_dir, experiment_name)
    output_prefix = output_prefix + "_" if not output_prefix.endswith("_") else output_prefix

    # Check if required parameters are provided
    if not metadata_file:
        logging.error("metadata_file is required.")
        raise ValueError("metadata_file is required.")
    if not input_fastq_dir:
        logging.error("input_fastq_dir is required.")
        raise ValueError("input_fastq_dir is required.")
    if not mode:
        logging.error("mode is required.")
        raise ValueError("mode is required.")
    if not output_dir:
        logging.error("output_dir is required.")
        raise ValueError("output_dir is required.")

    ## Debug logging
    # logging.info(f"metadata file: {metadata_file}")
    # logging.info(f"input fastq dir: {input_fastq_dir}")
    # logging.info(f"mode: {mode}")
    # logging.info(f"output dir: {output_dir}")
    # with open(metadata_file, 'r') as file:
    #     print(file.read())

    sample_names = pm.process_metadata(metadata_file, input_fastq_dir)[1]

    count_file_paths = [os.path.join(output_dir, f'{sample}_{mode}/{sample}_{mode}.count.txt') for sample in sample_names if sample and mode]

    logging.info("Building counts matrix...")
    sgrna_count_matrix = count_matrix(count_file_paths)

    save_matrix_with_excel(sgrna_count_matrix, output_prefix, "count")

    logging.info("Building CPM matrix...")
    sgrna_cpm_matrix = cpm_matrix(sgrna_count_matrix)
    save_matrix_with_excel(sgrna_cpm_matrix, output_prefix, "cpm")

    logging.info("Building CPM files for each sample...")
    for sample in sample_names:
        sample_cpm_df = sgrna_cpm_matrix[[sample]].copy()
        sample_cpm_df.columns = ['cpm']
        sample_cpm_df.to_csv(f"{output_dir}/{sample}_{mode}/{sample}_{mode}_cpm.txt", sep='\t')
        sample_cpm_df.to_csv(f"{output_dir}/{sample}_{mode}/{sample}_{mode}_cpm.csv")

    logging.info("Done.")

if __name__ == "__main__":

    logging.basicConfig(
        level=logging.INFO,
        format="[{asctime}] [{levelname}] {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler("output/cpm_matrix.log")
        ]
    )

    parser = argparse.ArgumentParser(description='Generate count and normalized count matrices.')
    parser.add_argument('--config', help='Path to the config file (default: config.yaml if exists, use "none" to disable config)')
    parser.add_argument('--metadata', help='Path to the metadata file')
    parser.add_argument('--input-fastq-dir', help='Path to the directory containing the FASTQ files')
    parser.add_argument('--output-prefix', help='Prefix for output files')
    parser.add_argument('--mode', default="fastq", help='Method used for counting (See config file)')
    args = parser.parse_args()
    count_and_cpm_matrix(
        config_file=args.config,
        metadata_file=args.metadata,
        input_fastq_dir=args.input_fastq_dir,
        output_prefix=args.output_prefix,
        mode=args.mode
    )
