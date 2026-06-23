import logging
import os
import subprocess

from conftest import (
    TEST_DATA_DIR,
    compare_fastq_files_by_read_id,
    compare_fastq_files_by_header,
    compare_fastq_files_read_order_preserved,
)

# run tests in the repo base dir using:
# pytest -sv test/test_reorient_fastq_parallel_2.py


def test_reorient_fastq_parallel_2_vs_expected(tmp_path):
    """parallel_2 output matches the pre-computed expected baseline (order-agnostic)."""
    input_r1 = os.path.join(TEST_DATA_DIR, 'test_sampled_R1.fastq.gz')
    input_r2 = os.path.join(TEST_DATA_DIR, 'test_sampled_R2.fastq.gz')
    expected_r1 = os.path.join(TEST_DATA_DIR, 'expected_R1.reoriented.fastq.gz')
    expected_r2 = os.path.join(TEST_DATA_DIR, 'expected_R2.reoriented.fastq.gz')

    output_r1 = tmp_path / 'v2p_output_R1.reoriented.fastq.gz'
    output_r2 = tmp_path / 'v2p_output_R2.reoriented.fastq.gz'

    args = [
        '--sequences-r1', 'GCTATTTCTAGCTCTAAAAC,TCCCACTCCTTTCAAGA',
        '--sequences-r2', 'GTTTTAGAGCTAGAAATAGC,GAAAGGACGAAACACCG',
        '--cpus', '4',
        '--batch-size', '10000',
        '--in-fastq-r1', input_r1,
        '--in-fastq-r2', input_r2,
        '--out-fastq-r1', str(output_r1),
        '--out-fastq-r2', str(output_r2),
    ]

    logging.info("Running reorient_fastq_parallel_2.py...")
    subprocess.run(['python', 'reorient_fastq_parallel_2.py'] + args,
                   check=True, capture_output=True, text=True)

    ok_r1, msg_r1 = compare_fastq_files_by_read_id(str(output_r1), expected_r1)
    assert ok_r1, f"R1 differs from expected: {msg_r1}"
    logging.info(f"✓ R1: {msg_r1}")

    ok_r2, msg_r2 = compare_fastq_files_by_read_id(str(output_r2), expected_r2)
    assert ok_r2, f"R2 differs from expected: {msg_r2}"
    logging.info(f"✓ R2: {msg_r2}")


def test_reorient_fastq_parallel_2_output_matches_v1(tmp_path):
    """parallel_2 and v1 produce identical results on the same input (order-agnostic, exact header)."""
    input_r1 = os.path.join(TEST_DATA_DIR, 'test_sampled_R1.fastq.gz')
    input_r2 = os.path.join(TEST_DATA_DIR, 'test_sampled_R2.fastq.gz')

    v1_r1 = tmp_path / 'v1_R1.reoriented.fastq.gz'
    v1_r2 = tmp_path / 'v1_R2.reoriented.fastq.gz'
    v2p_r1 = tmp_path / 'v2p_R1.reoriented.fastq.gz'
    v2p_r2 = tmp_path / 'v2p_R2.reoriented.fastq.gz'

    shared_args = [
        '--sequences-r1', 'GCTATTTCTAGCTCTAAAAC,TCCCACTCCTTTCAAGA',
        '--sequences-r2', 'GTTTTAGAGCTAGAAATAGC,GAAAGGACGAAACACCG',
        '--batch-size', '10000',
        '--in-fastq-r1', input_r1,
        '--in-fastq-r2', input_r2,
    ]

    logging.info("Running v1...")
    subprocess.run(['python', 'reorient_fastq_parallel.py'] + shared_args +
                   ['--cpus', '4', '--out-fastq-r1', str(v1_r1), '--out-fastq-r2', str(v1_r2)],
                   check=True, capture_output=True, text=True)

    logging.info("Running parallel_2...")
    subprocess.run(['python', 'reorient_fastq_parallel_2.py'] + shared_args +
                   ['--cpus', '4', '--out-fastq-r1', str(v2p_r1), '--out-fastq-r2', str(v2p_r2)],
                   check=True, capture_output=True, text=True)

    ok_r1, msg_r1 = compare_fastq_files_by_header(str(v1_r1), str(v2p_r1))
    assert ok_r1, f"R1 mismatch v1 vs parallel_2: {msg_r1}"
    logging.info(f"✓ R1: {msg_r1}")

    ok_r2, msg_r2 = compare_fastq_files_by_header(str(v1_r2), str(v2p_r2))
    assert ok_r2, f"R2 mismatch v1 vs parallel_2: {msg_r2}"
    logging.info(f"✓ R2: {msg_r2}")


def test_reorient_fastq_parallel_2_preserves_input_read_order(tmp_path):
    """parallel_2 preserves interleaved read order from input (pool.imap is ordered)."""
    input_r1 = os.path.join(TEST_DATA_DIR, 'test_sampled_R1.fastq.gz')
    input_r2 = os.path.join(TEST_DATA_DIR, 'test_sampled_R2.fastq.gz')

    output_r1 = tmp_path / 'order_v2p_R1.fastq.gz'
    output_r2 = tmp_path / 'order_v2p_R2.fastq.gz'

    args = [
        'python', 'reorient_fastq_parallel_2.py',
        '--sequences-r1', 'GCTATTTCTAGCTCTAAAAC,TCCCACTCCTTTCAAGA',
        '--sequences-r2', 'GTTTTAGAGCTAGAAATAGC,GAAAGGACGAAACACCG',
        '--cpus', '4',
        '--batch-size', '10000',
        '--in-fastq-r1', input_r1,
        '--in-fastq-r2', input_r2,
        '--out-fastq-r1', str(output_r1),
        '--out-fastq-r2', str(output_r2),
    ]

    logging.info("Running parallel_2...")
    subprocess.run(args, check=True, capture_output=True, text=True)

    ok, msg = compare_fastq_files_read_order_preserved(
        input_r1, input_r2, str(output_r1), str(output_r2))
    assert ok, f"Record order not preserved: {msg}"
    logging.info(f"✓ {msg}")
