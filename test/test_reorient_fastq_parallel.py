import os
import gzip
import subprocess
import logging
import hashlib
from Bio import SeqIO

# run tests in the repo base dir using:
# pytest -sv -o log_cli=true [optional: <specific test file or function>]

# non-test functions

def compare_fastq_files_strict(file1, file2):
    """Compare two FASTQ files for exact identity (record-by-record)"""
    with gzip.open(file1, 'rt') as f1, gzip.open(file2, 'rt') as f2:
        records1 = list(SeqIO.parse(f1, "fastq"))
        records2 = list(SeqIO.parse(f2, "fastq"))

    if len(records1) != len(records2):
        return False, f"Different number of records: {len(records1)} vs {len(records2)}"

    for i, (r1, r2) in enumerate(zip(records1, records2)):
        if r1.id != r2.id:
            return False, f"Record {i}: ID mismatch: {r1.id} vs {r2.id}"
        if str(r1.seq) != str(r2.seq):
            return False, f"Record {i}: Sequence mismatch"
        if r1.letter_annotations["phred_quality"] != r2.letter_annotations["phred_quality"]:
            return False, f"Record {i}: Quality mismatch"

    return True, "Files are identical"


def _strip_read_suffix(read_id):
    """Remove /1 or /2 suffix from read ID (used in paired-end sequencing)."""
    if read_id.endswith('/1') or read_id.endswith('/2'):
        return read_id[:-2]
    return read_id


def compare_fastq_files_by_read_id(file1, file2):
    """Compare FASTQ files by base read ID (order-agnostic, header-lenient).

    Extracts base read ID by stripping /1 or /2 suffixes and matches records on that.
    Allows headers to differ in their suffix formatting while validating sequences/qualities match.

    Returns (True/False, message) tuple.
    """
    dict1 = {}
    dict2 = {}

    with gzip.open(file1, 'rt') as f:
        for record in SeqIO.parse(f, "fastq"):
            key = _strip_read_suffix(record.id)
            dict1[key] = record

    with gzip.open(file2, 'rt') as f:
        for record in SeqIO.parse(f, "fastq"):
            key = _strip_read_suffix(record.id)
            dict2[key] = record

    if len(dict1) != len(dict2):
        return False, f"Different number of unique base read IDs: {len(dict1)} vs {len(dict2)}"

    if set(dict1.keys()) != set(dict2.keys()):
        only_in_1 = set(dict1.keys()) - set(dict2.keys())
        only_in_2 = set(dict2.keys()) - set(dict1.keys())
        msg = f"Different base read IDs"
        if only_in_1:
            msg += f"\nOnly in file1: {list(only_in_1)[:3]}..."
        if only_in_2:
            msg += f"\nOnly in file2: {list(only_in_2)[:3]}..."
        return False, msg

    for base_id in dict1:
        r1 = dict1[base_id]
        r2 = dict2[base_id]
        if str(r1.seq) != str(r2.seq):
            return False, f"Record {base_id}: Sequence mismatch"
        if r1.letter_annotations["phred_quality"] != r2.letter_annotations["phred_quality"]:
            return False, f"Record {base_id}: Quality mismatch"

    return True, f"All {len(dict1)} records match (order-agnostic, base read ID)"


def compare_fastq_files_by_header(file1, file2):
    """Compare FASTQ files by exact header match (order-agnostic, header-strict).

    Requires exact full header ID match (including /1 or /2 suffixes) along with
    sequences and qualities. This validates not just reorientation logic but also
    that headers are formatted identically.

    Returns (True/False, message) tuple.
    """
    dict1 = {}
    dict2 = {}

    with gzip.open(file1, 'rt') as f:
        for record in SeqIO.parse(f, "fastq"):
            dict1[record.id] = record

    with gzip.open(file2, 'rt') as f:
        for record in SeqIO.parse(f, "fastq"):
            dict2[record.id] = record

    if len(dict1) != len(dict2):
        return False, f"Different number of unique headers: {len(dict1)} vs {len(dict2)}"

    if set(dict1.keys()) != set(dict2.keys()):
        only_in_1 = set(dict1.keys()) - set(dict2.keys())
        only_in_2 = set(dict2.keys()) - set(dict1.keys())
        msg = f"Different headers"
        if only_in_1:
            msg += f"\nOnly in file1: {list(only_in_1)[:3]}..."
        if only_in_2:
            msg += f"\nOnly in file2: {list(only_in_2)[:3]}..."
        return False, msg

    for header in dict1:
        r1 = dict1[header]
        r2 = dict2[header]
        if str(r1.seq) != str(r2.seq):
            return False, f"Record {header}: Sequence mismatch"
        if r1.letter_annotations["phred_quality"] != r2.letter_annotations["phred_quality"]:
            return False, f"Record {header}: Quality mismatch"

    return True, f"All {len(dict1)} records match (order-agnostic, exact header)"


def compare_fastq_files_read_order_preserved(input_r1, input_r2, output_r1, output_r2):
    """Compare that output preserves the order of base read IDs from input.

    Creates interleaved sequences of R1 and R2 (R1[0], R2[0], R1[1], R2[1], ...)
    for both input and output, then verifies that the base read IDs appear in the
    same order. This validates that reorientation preserves record order without
    reordering them.

    Only checks that base read IDs match at each position; doesn't validate
    sequence or quality content.

    Returns (True/False, message) tuple.
    """
    input_ids = []
    with gzip.open(input_r1, 'rt') as f1, gzip.open(input_r2, 'rt') as f2:
        records_r1 = list(SeqIO.parse(f1, "fastq"))
        records_r2 = list(SeqIO.parse(f2, "fastq"))

    if len(records_r1) != len(records_r2):
        return False, f"Input files have different record counts: R1={len(records_r1)}, R2={len(records_r2)}"

    for r1, r2 in zip(records_r1, records_r2):
        input_ids.append(_strip_read_suffix(r1.id))
        input_ids.append(_strip_read_suffix(r2.id))

    output_ids = []
    with gzip.open(output_r1, 'rt') as f1, gzip.open(output_r2, 'rt') as f2:
        out_records_r1 = list(SeqIO.parse(f1, "fastq"))
        out_records_r2 = list(SeqIO.parse(f2, "fastq"))

    if len(out_records_r1) != len(out_records_r2):
        return False, f"Output files have different record counts: R1={len(out_records_r1)}, R2={len(out_records_r2)}"

    for r1, r2 in zip(out_records_r1, out_records_r2):
        output_ids.append(_strip_read_suffix(r1.id))
        output_ids.append(_strip_read_suffix(r2.id))

    if len(input_ids) != len(output_ids):
        return False, f"Different total record counts: input={len(input_ids)}, output={len(output_ids)}"

    for i, (input_id, output_id) in enumerate(zip(input_ids, output_ids)):
        if input_id != output_id:
            return False, f"Position {i}: Base read ID mismatch: {input_id} (input) vs {output_id} (output)"

    return True, f"All {len(input_ids)} records match in order (order-sensitive, base read ID only)"


# test functions

def test_reorient_fastq_parallel_v1_vs_expected(tmp_path):
    """
    Test that reorient_fastq_parallel.py (v1) produces expected output.
    This validates against the baseline expected output files:
    - expected_R1.reoriented.fastq.gz
    - expected_R2.reoriented.fastq.gz
    """
    logging.basicConfig(level=logging.INFO,
                        format='[%(asctime)s] [%(levelname)s] %(message)s',
                        handlers=[logging.StreamHandler()])

    logging.info("Running test_reorient_fastq_parallel_v1_vs_expected")

    test_data_dir = os.path.join(os.path.dirname(__file__), 'test_data')
    input_r1 = os.path.join(test_data_dir, 'test_sampled_R1.fastq.gz')
    input_r2 = os.path.join(test_data_dir, 'test_sampled_R2.fastq.gz')
    expected_r1 = os.path.join(test_data_dir, 'expected_R1.reoriented.fastq.gz')
    expected_r2 = os.path.join(test_data_dir, 'expected_R2.reoriented.fastq.gz')

    output_r1 = tmp_path / 'v1_output_R1.reoriented.fastq.gz'
    output_r2 = tmp_path / 'v1_output_R2.reoriented.fastq.gz'

    print(f"Temporary path: {tmp_path}", flush=True)

    # IMPORTANT: use the same args as were used to create expected files!
    # changing --cpus will cause strict comparison to fail because of record ordering,
    # but a lenient comparison should still pass.
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

    logging.info(f"Running reorient_fastq_parallel.py (v1) with args: {args}")
    result = subprocess.run(['python', 'reorient_fastq_parallel.py'] + args,
                          check=True, capture_output=True, text=True)
    logging.info(f"Script finished with return code: {result.returncode}")

    logging.info("Comparing output with expected files...")
    with gzip.open(output_r1, 'rt') as f1, gzip.open(expected_r1, 'rt') as f2:
        output_r1_content = f1.read()
        expected_r1_content = f2.read()
        if output_r1_content != expected_r1_content:
            print("\n" + "="*70)
            print("R1 OUTPUT DIFFERS FROM EXPECTED - Showing first 3 headers from each:")
            print("="*70)
            with gzip.open(output_r1, 'rt') as f:
                print("\nACTUAL (from v1 output):")
                for i, line in enumerate(f):
                    if i < 12:
                        print(line.rstrip())
                    else:
                        break
            with gzip.open(expected_r1, 'rt') as f:
                print("\nEXPECTED (from baseline):")
                for i, line in enumerate(f):
                    if i < 12:
                        print(line.rstrip())
                    else:
                        break
            print("="*70)
        assert output_r1_content == expected_r1_content, \
            "R1 output file content does not match expected content"
    logging.info("R1 matches expected output.")

    with gzip.open(output_r2, 'rt') as f1, gzip.open(expected_r2, 'rt') as f2:
        output_r2_content = f1.read()
        expected_r2_content = f2.read()
        if output_r2_content != expected_r2_content:
            print("\n" + "="*70)
            print("R2 OUTPUT DIFFERS FROM EXPECTED - Showing first 3 headers from each:")
            print("="*70)
            with gzip.open(output_r2, 'rt') as f:
                print("\nACTUAL (from v1 output):")
                for i, line in enumerate(f):
                    if i < 12:
                        print(line.rstrip())
                    else:
                        break
            with gzip.open(expected_r2, 'rt') as f:
                print("\nEXPECTED (from baseline):")
                for i, line in enumerate(f):
                    if i < 12:
                        print(line.rstrip())
                    else:
                        break
            print("="*70)
        assert output_r2_content == expected_r2_content, \
            "R2 output file content does not match expected content"
    logging.info("R2 matches expected output.")

    logging.info("="*60)
    logging.info("SUCCESS: v1 output matches expected baseline!")
    logging.info("="*60)
