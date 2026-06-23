import gzip
import logging
import os
import subprocess

from conftest import TEST_DATA_DIR

# run tests in the repo base dir using:
# pytest -sv test/test_reorient_fastq_parallel.py
#
# NOTE: reorient_fastq_parallel.py is the original v1 script (archival).
# This test only needs to run if v1 code is changed, which is not expected.


def test_reorient_fastq_parallel_v1_vs_expected(tmp_path):
    """v1 produces byte-identical output to the pre-computed expected baseline."""
    input_r1 = os.path.join(TEST_DATA_DIR, 'test_sampled_R1.fastq.gz')
    input_r2 = os.path.join(TEST_DATA_DIR, 'test_sampled_R2.fastq.gz')
    expected_r1 = os.path.join(TEST_DATA_DIR, 'expected_R1.reoriented.fastq.gz')
    expected_r2 = os.path.join(TEST_DATA_DIR, 'expected_R2.reoriented.fastq.gz')

    output_r1 = tmp_path / 'v1_output_R1.reoriented.fastq.gz'
    output_r2 = tmp_path / 'v1_output_R2.reoriented.fastq.gz'

    # IMPORTANT: --cpus must match what was used to generate the expected files,
    # otherwise record ordering will differ (parallel batches may reorder).
    args = [
        '--sequences-r1', 'GCTATTTCTAGCTCTAAAAC,TCCCACTCCTTTCAAGA',
        '--sequences-r2', 'GTTTTAGAGCTAGAAATAGC,GAAAGGACGAAACACCG',
        '--cpus', '16',
        '--batch-size', '10000',
        '--in-fastq-r1', input_r1,
        '--in-fastq-r2', input_r2,
        '--out-fastq-r1', str(output_r1),
        '--out-fastq-r2', str(output_r2),
    ]

    logging.info("Running reorient_fastq_parallel.py (v1)...")
    subprocess.run(['python', 'reorient_fastq_parallel.py'] + args,
                   check=True, capture_output=True, text=True)

    for out_path, exp_path, label in [
        (output_r1, expected_r1, 'R1'),
        (output_r2, expected_r2, 'R2'),
    ]:
        with gzip.open(out_path, 'rt') as f_out, gzip.open(exp_path, 'rt') as f_exp:
            out_content = f_out.read()
            exp_content = f_exp.read()
        if out_content != exp_content:
            with gzip.open(out_path, 'rt') as f:
                actual_head = [next(f) for _ in range(12)]
            with gzip.open(exp_path, 'rt') as f:
                expected_head = [next(f) for _ in range(12)]
            logging.error(f"{label} actual (first 3 records):\n{''.join(actual_head)}")
            logging.error(f"{label} expected (first 3 records):\n{''.join(expected_head)}")
        assert out_content == exp_content, f"{label} output does not match expected"
        logging.info(f"✓ {label} matches expected baseline")
