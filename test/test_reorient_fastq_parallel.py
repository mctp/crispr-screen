import os
import gzip
import subprocess
import logging

# run tests in the repo base dir
# pytest -sv -o log_cli=true (OPTIONAL: add specific test script to end of command)

def test_reorient_fastq_parallel(tmp_path):
    # Configure logging to output to the console
    logging.basicConfig(level=logging.INFO, 
                        format='[%(asctime)s] [%(levelname)s] %(message)s',
                        handlers=[logging.StreamHandler()])
    
    logging.info("Running test_reorient_fastq_parallel")

    # Define paths to test data
    test_data_dir = os.path.join(os.path.dirname(__file__), 'test_data')
    input_r1 = os.path.join(test_data_dir, 'test_sampled_R1.fastq.gz')
    input_r2 = os.path.join(test_data_dir, 'test_sampled_R2.fastq.gz')
    expected_r1 = os.path.join(test_data_dir, 'expected_R1.reoriented.fastq.gz')
    expected_r2 = os.path.join(test_data_dir, 'expected_R2.reoriented.fastq.gz')
    # Define paths to output data
    output_r1 = tmp_path / 'output_R1.reoriented.fastq.gz'
    output_r2 = tmp_path / 'output_R2.reoriented.fastq.gz'
    # interleaved_gz = tmp_path / 'interleaved.fastq.gz'
    plot_prefix = tmp_path / 'plot'

    # Print the temporary path
    print(f"Temporary path: {tmp_path}", flush=True)

    # Run the main function with test arguments
    args = [
        '--sequences-r1', 'GCTATTTCTAGCTCTAAAAC,TCCCACTCCTTTCAAGA',
        '--sequences-r2', 'GTTTTAGAGCTAGAAATAGC,GAAAGGACGAAACACCG',
        '--cpus', '16',
        '--batch-size', '1000',
        '--in-fastq-r1', input_r1,
        '--in-fastq-r2', input_r2,
        '--out-fastq-r1', str(output_r1),
        '--out-fastq-r2', str(output_r2),
        '--validate',
        '--plot',
        '--plot-prefix', str(plot_prefix)
    ]
    
    logging.info(f"Running subprocess with args: {args}")
    result = subprocess.run(['python', 'reorient_fastq_parallel.py'] + args, check=True)
    logging.info(f"Subprocess finished with return code: {result.returncode}")

    #  Compare the uncompressed contents of the output files with the expected files
    with gzip.open(output_r1, 'rt') as f1, gzip.open(expected_r1, 'rt') as f2:
        output_r1_content = f1.read()
        expected_r1_content = f2.read()
        assert output_r1_content == expected_r1_content, "R1 output file content does not match expected content"

    with gzip.open(output_r2, 'rt') as f1, gzip.open(expected_r2, 'rt') as f2:
        output_r2_content = f1.read()
        expected_r2_content = f2.read()
        assert output_r2_content == expected_r2_content, "R2 output file content does not match expected content"
