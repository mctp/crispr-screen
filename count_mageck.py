import os
import subprocess
import argparse

def run_mageck_count(mode, sample, output_dir, references_dir, sgrna_list_name, in_r1_fastq, in_r2_fastq, trim5_lengths="", keep_tmp=False, save_unmapped=False, search_revcomp=False, ignore_r2=False):
    sample_dir = os.path.join(output_dir, f"{sample}_{mode}" if mode else sample)
    sample_prefix = os.path.join(sample_dir, f"{sample}_{mode}" if mode else sample)
    out_bam = f"{sample_prefix}.aligned.bam"
    mageck_out_prefix = sample_prefix

    os.makedirs(sample_dir, exist_ok=True)

    if mode in ["bam", ""]:
        sgrna_list_prefix = os.path.join(references_dir, sgrna_list_name)
        sgrna_list_file = f"{sgrna_list_prefix}.txt"
        os.makedirs(references_dir, exist_ok=True)
        subprocess.run([
            "mageck", "count",
            "-l", sgrna_list_file,
            "-n", mageck_out_prefix,
            "--sample-label", sample,
            "--fastq", out_bam
        ])
    elif mode == "fastq":
        sgrna_list_file = f"{sgrna_list_name}.txt"
        other_opts = []

        if trim5_lengths:
            other_opts.extend(["--trim-5", trim5_lengths])
        if keep_tmp:
            other_opts.append("--keep-tmp")
        if save_unmapped:
            other_opts.append("--unmapped-to-file")
        if trim5_lengths or keep_tmp or save_unmapped:
            print("Forcing sgrna length to 20 because options trim5_lengths, keep_tmp, and/or save_unmapped are set in config.")
            other_opts.extend(["--sgrna-len", "20"])
        if not ignore_r2:
            other_opts.extend(["--fastq-2", in_r2_fastq, "--count-pair", "True"])
        if search_revcomp:
            other_opts.append("--reverse-complement")

        subprocess.run([
            "mageck", "count",
            "-l", sgrna_list_file,
            "-n", mageck_out_prefix,
            "--sample-label", sample,
            "--fastq", in_r1_fastq
        ] + other_opts)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Run MAGeCK count with specified parameters.")
    parser.add_argument("--mode", required=True, help="Mode to run MAGeCK count (bam or fastq).")
    parser.add_argument("--sample", required=True, help="Sample name.")
    parser.add_argument("--output_dir", required=True, help="Output directory.")
    parser.add_argument("--references_dir", required=True, help="References directory.")
    parser.add_argument("--sgrna_list_name", required=True, help="sgRNA list name.")
    parser.add_argument("--in_r1_fastq", required=True, help="Input R1 FASTQ file.")
    parser.add_argument("--in_r2_fastq", required=False, help="Input R2 FASTQ file.")
    parser.add_argument("--trim5_lengths", default="", help="Trim 5 lengths.")
    parser.add_argument("--keep_tmp", action="store_true", help="Keep temporary files.")
    parser.add_argument("--save_unmapped", action="store_true", help="Save unmapped reads.")
    parser.add_argument("--search_revcomp", action="store_true", help="Search reverse complement.")
    parser.add_argument("--ignore_r2", action="store_true", help="Ignore R2 FASTQ file.")

    args = parser.parse_args()

    run_mageck_count(
        mode=args.mode,
        sample=args.sample,
        output_dir=args.output_dir,
        references_dir=args.references_dir,
        sgrna_list_name=args.sgrna_list_name,
        in_r1_fastq=args.in_r1_fastq,
        in_r2_fastq=args.in_r2_fastq,
        trim5_lengths=args.trim5_lengths,
        keep_tmp=args.keep_tmp,
        save_unmapped=args.save_unmapped,
        search_revcomp=args.search_revcomp,
        ignore_r2=args.ignore_r2
    )
