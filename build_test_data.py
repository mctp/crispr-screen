import argparse
import subprocess
import sys
import gzip
import shutil

def check_seqtk():
    try:
        result = subprocess.run(["seqtk"], capture_output=True, text=True)
        output = result.stdout + result.stderr
        if "Version" not in output:
            raise FileNotFoundError
        version_line = [line for line in output.split('\n') if "Version" in line][0]
        version = version_line.split()[1]
        print(f"seqtk version: {version}")
    except (FileNotFoundError, IndexError):
        print("Error: seqtk is not available or the version could not be determined. Please install seqtk and try again.")
        sys.exit(1)

def check_pigz():
    try:
        result = subprocess.run(["pigz", "--version"], capture_output=True, text=True)
        if result.returncode == 0:
            version_line = result.stdout.split('\n')[0]
            version = version_line.split()[1]
            print(f"pigz version: {version}")
            return True
        else:
            raise FileNotFoundError
    except FileNotFoundError:
        print("pigz is not available. Falling back to gzip.")
        return False

def combine_fastq_files(fastq_files, output_file):
    cmd = f"cat {' '.join(fastq_files)} > {output_file}"
    subprocess.run(cmd, shell=True, check=True)

def sample_fastq(input_file, output_file, nreads, seed, use_pigz):
    temp_output_file = output_file.replace('.gz', '')
    cmd = f"seqtk sample -s{seed} {input_file} {nreads} > {temp_output_file}"
    subprocess.run(cmd, shell=True, check=True)
    if use_pigz:
        cmd = f"pigz -c {temp_output_file} > {output_file}"
    else:
        with open(temp_output_file, 'rb') as f_in, gzip.open(output_file, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    subprocess.run(cmd, shell=True, check=True)
    subprocess.run(f"rm {temp_output_file}", shell=True, check=True)

def is_gzipped(file):
    with open(file, 'rb') as f:
        return f.read(2) == b'\x1f\x8b'

def main():
    check_seqtk()
    use_pigz = check_pigz()
    
    parser = argparse.ArgumentParser(description="Combine and sample paired-end fastq files.")
    parser.add_argument("-1", "--fastq1", required=True, help="Comma-delimited list of R1 fastq files.")
    parser.add_argument("-2", "--fastq2", required=True, help="Comma-delimited list of R2 fastq files.")
    parser.add_argument("-n", "--nreads", type=int, default=10000, help="Number of reads to sample.")
    parser.add_argument("-s", "--seed", type=int, default=1234, help="Seed for random sampling.")
    parser.add_argument("-o", "--output", required=True, help="Output prefix for combined and sampled files.")
    
    args = parser.parse_args()
    
    fastq1_files = args.fastq1.split(',')
    fastq2_files = args.fastq2.split(',')
    
    if all(is_gzipped(f) for f in fastq1_files):
        combined_fastq1 = f"{args.output}_combined_R1.fastq.gz"
    else:
        combined_fastq1 = f"{args.output}_combined_R1.fastq"

    if all(is_gzipped(f) for f in fastq2_files):
        combined_fastq2 = f"{args.output}_combined_R2.fastq.gz"
    else:
        combined_fastq2 = f"{args.output}_combined_R2.fastq"
    print("Combining fastq files...")
    combine_fastq_files(fastq1_files, combined_fastq1)
    combine_fastq_files(fastq2_files, combined_fastq2)
    
    sampled_fastq1 = f"{args.output}_sampled_R1.fastq.gz"
    sampled_fastq2 = f"{args.output}_sampled_R2.fastq.gz"
    
    print(f"Sampling fastq files (n = {args.nreads}, seed = {args.seed})...")
    sample_fastq(combined_fastq1, sampled_fastq1, args.nreads, args.seed, use_pigz)
    sample_fastq(combined_fastq2, sampled_fastq2, args.nreads, args.seed, use_pigz)

if __name__ == "__main__":
    main()
