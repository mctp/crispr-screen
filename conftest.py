import gzip
import logging
import os

import pytest
from Bio import SeqIO

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'test', 'test_data')


@pytest.fixture(autouse=True, scope='session')
def configure_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='[%(asctime)s] [%(levelname)s] %(message)s',
        handlers=[logging.StreamHandler()]
    )


def _strip_read_suffix(read_id):
    if read_id.endswith('/1') or read_id.endswith('/2'):
        return read_id[:-2]
    return read_id


def compare_fastq_files_strict(file1, file2):
    """Record-by-record comparison (order-sensitive, exact match)."""
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


def compare_fastq_files_by_read_id(file1, file2):
    """Order-agnostic comparison by base read ID (strips /1 /2 suffixes)."""
    dict1 = {}
    dict2 = {}

    with gzip.open(file1, 'rt') as f:
        for record in SeqIO.parse(f, "fastq"):
            dict1[_strip_read_suffix(record.id)] = record

    with gzip.open(file2, 'rt') as f:
        for record in SeqIO.parse(f, "fastq"):
            dict2[_strip_read_suffix(record.id)] = record

    if len(dict1) != len(dict2):
        return False, f"Different number of unique base read IDs: {len(dict1)} vs {len(dict2)}"

    if set(dict1.keys()) != set(dict2.keys()):
        only_in_1 = set(dict1.keys()) - set(dict2.keys())
        only_in_2 = set(dict2.keys()) - set(dict1.keys())
        msg = "Different base read IDs"
        if only_in_1:
            msg += f"\n  Only in file1: {list(only_in_1)[:3]}..."
        if only_in_2:
            msg += f"\n  Only in file2: {list(only_in_2)[:3]}..."
        return False, msg

    for base_id, r1 in dict1.items():
        r2 = dict2[base_id]
        if str(r1.seq) != str(r2.seq):
            return False, f"Record {base_id}: Sequence mismatch"
        if r1.letter_annotations["phred_quality"] != r2.letter_annotations["phred_quality"]:
            return False, f"Record {base_id}: Quality mismatch"

    return True, f"All {len(dict1)} records match (order-agnostic, base read ID)"


def compare_fastq_files_by_header(file1, file2):
    """Order-agnostic comparison by exact full header (no suffix stripping)."""
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
        msg = "Different headers"
        if only_in_1:
            msg += f"\n  Only in file1: {list(only_in_1)[:3]}..."
        if only_in_2:
            msg += f"\n  Only in file2: {list(only_in_2)[:3]}..."
        return False, msg

    for header, r1 in dict1.items():
        r2 = dict2[header]
        if str(r1.seq) != str(r2.seq):
            return False, f"Record {header}: Sequence mismatch"
        if r1.letter_annotations["phred_quality"] != r2.letter_annotations["phred_quality"]:
            return False, f"Record {header}: Quality mismatch"

    return True, f"All {len(dict1)} records match (order-agnostic, exact header)"


def compare_fastq_files_read_order_preserved(input_r1, input_r2, output_r1, output_r2):
    """Verify interleaved read order (R1[0],R2[0],R1[1],...) is preserved from input to output."""
    with gzip.open(input_r1, 'rt') as f1, gzip.open(input_r2, 'rt') as f2:
        in_r1 = list(SeqIO.parse(f1, "fastq"))
        in_r2 = list(SeqIO.parse(f2, "fastq"))

    if len(in_r1) != len(in_r2):
        return False, f"Input files have different record counts: R1={len(in_r1)}, R2={len(in_r2)}"

    input_ids = []
    for r1, r2 in zip(in_r1, in_r2):
        input_ids.append(_strip_read_suffix(r1.id))
        input_ids.append(_strip_read_suffix(r2.id))

    with gzip.open(output_r1, 'rt') as f1, gzip.open(output_r2, 'rt') as f2:
        out_r1 = list(SeqIO.parse(f1, "fastq"))
        out_r2 = list(SeqIO.parse(f2, "fastq"))

    if len(out_r1) != len(out_r2):
        return False, f"Output files have different record counts: R1={len(out_r1)}, R2={len(out_r2)}"

    output_ids = []
    for r1, r2 in zip(out_r1, out_r2):
        output_ids.append(_strip_read_suffix(r1.id))
        output_ids.append(_strip_read_suffix(r2.id))

    if len(input_ids) != len(output_ids):
        return False, f"Different total record counts: input={len(input_ids)}, output={len(output_ids)}"

    for i, (in_id, out_id) in enumerate(zip(input_ids, output_ids)):
        if in_id != out_id:
            return False, f"Position {i}: {in_id} (input) vs {out_id} (output)"

    return True, f"All {len(input_ids)} records match in order"
