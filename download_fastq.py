#!/usr/bin/env python3
"""
Download FASTQ files listed in sample_metadata.yaml from Google Cloud Storage.
Parses the flowcell from each filename and constructs the GCS path:
  gs://<bucket>/<flowcell>/<filename>
If the flowcell directory is not found, falls back through alt_fastq_dirs (from config).
Modify FLOWCELL_PATTERN to match your actual filename convention.
"""

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys

import yaml

GCLOUD_FALLBACK_PATHS = [
    '/home/mambauser/google-cloud-sdk/bin/gcloud',
    '/root/google-cloud-sdk/bin/gcloud',
    '/usr/lib/google-cloud-sdk/bin/gcloud',
    '/snap/bin/gcloud',
]

FLOWCELL_PATTERN = re.compile(r'^mctp_SI_\d+_([^_]+)_\d+_\d+\.fq\.gz$')


def resolve_gcloud(cmd_string):
    """Replace 'gcloud' with its full path when it is not on PATH (e.g. inside Docker)."""
    parts = cmd_string.split()
    if parts[0] != 'gcloud':
        return parts
    gcloud_bin = shutil.which('gcloud')
    if not gcloud_bin:
        for path in GCLOUD_FALLBACK_PATHS:
            if os.path.isfile(path) and os.access(path, os.X_OK):
                gcloud_bin = path
                break
    if not gcloud_bin:
        logging.error("gcloud not found. Add it to PATH or install the Google Cloud SDK.")
        sys.exit(1)
    return [gcloud_bin] + parts[1:]


def parse_args():
    parser = argparse.ArgumentParser(
        description='Download FASTQ files from GCS based on config.yaml and sample_metadata.yaml',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--config', default='config.yaml', help='Pipeline config file')
    parser.add_argument('--dry-run', action='store_true', help='Print commands without executing')
    parser.add_argument('--skip-file-check', action='store_true',
                        help='Skip verifying file existence in the bucket before downloading')
    parser.add_argument('--project', default='',
                        help='GCP project ID to verify against the active gcloud project; '
                             'overrides gcloud_project in config')
    return parser.parse_args()


def load_yaml(path):
    with open(path) as f:
        return yaml.safe_load(f)


def collect_fastq_files(metadata_file):
    """Return a flat list of all FASTQ file entries from the metadata."""
    data = load_yaml(metadata_file)
    if 'samples' not in data:
        raise ValueError(f"'samples' key not found in {metadata_file}")
    files = []
    for sample in data['samples']:
        if not isinstance(sample, dict):
            continue
        if sample.get('sample', '').startswith('#'):
            continue
        for key in ('fastq_r1', 'fastq_r2'):
            entries = sample.get(key, [])
            if isinstance(entries, str):
                entries = [entries]
            files.extend(entries)
    return files


def flowcell_from_filename(filename):
    m = FLOWCELL_PATTERN.match(os.path.basename(filename))
    if not m:
        raise ValueError(f"Cannot parse flowcell from filename: {filename}")
    return m.group(1)


def gcs_file_exists(ls_cmd, path):
    result = subprocess.run(resolve_gcloud(ls_cmd) + [path], capture_output=True)
    return result.returncode == 0


def find_gcs_path(ls_cmd, gcs_base, flowcell, alt_dirs, fname):
    """Try flowcell subdir first, then each alt_fastq_dirs entry."""
    candidates = [f"{gcs_base}/{flowcell}/{fname}"] + [f"{gcs_base}/{d}/{fname}" for d in alt_dirs]
    for path in candidates:
        logging.info(f"  checking {path} ...")
        if gcs_file_exists(ls_cmd, path):
            return path
    return None


def validate_gcp_project(expected_project):
    """Exit with an error if the active gcloud project does not match expected_project."""
    active_project = (
        os.environ.get("GCLOUD_PROJECT")
        or os.environ.get("CLOUDSDK_CORE_PROJECT")
        or ""
    )
    if not active_project:
        result = subprocess.run(
            resolve_gcloud("gcloud config get-value project"),
            capture_output=True, text=True,
        )
        lines = [l.strip() for l in (result.stdout + result.stderr).splitlines()
                 if l.strip() and not l.startswith("Your active configuration")]
        active_project = lines[-1] if lines else ""
    if active_project != expected_project:
        logging.error(
            f"Active gcloud project '{active_project}' does not match "
            f"expected '{expected_project}'.\n"
            f"  To switch: gcloud config set project {expected_project}"
        )
        sys.exit(1)
    logging.info(f"GCP project verified: {active_project}")


def main():
    logging.basicConfig(
        level=logging.INFO,
        format="[{asctime}] [{levelname}] {message}",
        style="{",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    args = parse_args()
    config = load_yaml(args.config)

    expected_project = args.project or config.get('gcloud_project', '') or ''
    if expected_project:
        validate_gcp_project(expected_project)
    else:
        logging.warning(
            "gcloud_project is not set in config.yaml and --project was not provided. "
            "Ensure the correct project is active: gcloud config set project <project-id>"
        )

    use_docker = config.get('docker_paths', False)
    dest_dir = config.get('docker_fastq_dir' if use_docker else 'fastq_dir', 'input')
    os.makedirs(dest_dir, exist_ok=True)

    metadata_file = config.get('metadata_file', 'sample_metadata.yaml')
    gcs_base = config.get('cloud_fastq_base_path', '').rstrip('/*').rstrip('/')
    ls_cmd = config.get('cloud_storage_ls_command', 'gcloud storage ls')
    cp_cmd = config.get('cloud_storage_cp_command', 'gcloud storage cp')
    alt_dirs = config.get('alt_fastq_dirs') or []

    if not gcs_base:
        logging.error("cloud_fastq_base_path is not set in config.yaml.")
        sys.exit(1)

    fastq_files = collect_fastq_files(metadata_file)
    if not fastq_files:
        logging.error("No FASTQ files found in metadata.")
        sys.exit(1)

    logging.info(f"Found {len(fastq_files)} FASTQ file(s) -> {dest_dir}/")
    if alt_dirs:
        logging.info(f"Fallback dirs: {alt_dirs}")

    failed = []
    for fpath in fastq_files:
        fname = os.path.basename(fpath)
        dst = os.path.join(dest_dir, fname)
        if os.path.exists(dst):
            logging.info(f"  already exists, skipping: {fname}")
            continue

        if args.skip_file_check:
            try:
                flowcell = flowcell_from_filename(fname)
                src = f"{gcs_base}/{flowcell}/{fname}"
            except ValueError:
                src = f"{gcs_base}/{fname}"
        else:
            try:
                flowcell = flowcell_from_filename(fname)
            except ValueError as e:
                logging.error(f"  {e}")
                failed.append(fname)
                continue
            src = find_gcs_path(ls_cmd, gcs_base, flowcell, alt_dirs, fname)
            if src is None:
                logging.error(f"  NOT FOUND: {fname} (tried flowcell={flowcell}, alt_dirs={alt_dirs})")
                failed.append(fname)
                continue

        logging.info(f"  found: {src}")
        cmd = resolve_gcloud(cp_cmd) + [src, dest_dir + '/']
        if args.dry_run:
            logging.info(f"  [dry-run] {' '.join(cmd)}")
        else:
            subprocess.run(cmd, check=True)

    if failed:
        logging.error(f"{len(failed)} file(s) not found: {failed}")
        sys.exit(1)

    logging.info("Done.")


if __name__ == '__main__':
    main()
