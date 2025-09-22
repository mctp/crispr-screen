import os
import sys
import argparse
import logging
import time
# sourced python scripts
import utilities.utilities as util
import utilities.process_metadata as pm
import reorient_fastq_parallel as reorient
import count_mageck as count
import count_mageck2 as count2
import count_reco as count_reco
import count_grep as cg
import cpm_matrix as cpm
import find_top_sequences as top
import distance_to_pattern_frequencies as dist
import mageck_analysis as ma

logging.basicConfig(
    level=logging.INFO,
    format="[{asctime}] [{levelname}] {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[
        logging.StreamHandler(sys.stdout),
        logging.FileHandler("output/run_pipeline.log")
    ]
)

def main():
    
    config = util.load_config(quiet=True)

    # logging.info("Loaded configuration:")
    # for key, value in config.items():
    #     logging.info(f"{key}: {value}")

    # establish defaults where config file is incomplete
    mode = config.get("mode", config.get("MODE", "fastq"))
    count_method_config = config.get("count_method", "mageck")
    
    # Handle multiple counting methods
    if isinstance(count_method_config, list):
        count_methods = count_method_config
    elif isinstance(count_method_config, str):
        count_methods = [count_method_config]
    else:
        logging.error(f"Invalid count_method configuration: {count_method_config}")
        sys.exit(1)
    
    # Validate all count methods
    valid_methods = ["mageck", "reco"]
    for method in count_methods:
        if method not in valid_methods:
            logging.error(f"Unknown count_method: {method}. Supported methods: {valid_methods}")
            sys.exit(1)

    search_revcomp = config.get("search_revcomp", config.get("SEARCH_REVCOMP", "FALSE"))

    reformat_sgrna_list = config.get("reformat_sgrna_list", "FALSE")
    skip_references = config.get("skip_references", "FALSE")
    force_references = config.get("force_references", "FALSE")
    skip_fastq_processing = config.get("skip_fastq_processing", "FALSE")
    skip_matrix = config.get("skip_matrix", "FALSE")
    skip_qc = config.get("skip_qc", "FALSE")
    skip_analysis = config.get("skip_analysis", "FALSE")
    analysis_only = config.get("analysis_only", "FALSE")
    paired_end_fastq = config.get("paired_end_fastq", "TRUE")
    reorient_fastq = config.get("reorient_fastq", "FALSE")
    sgrna_list_name = config.get("sgrna_list_name", config.get("SGRNA_LIST_NAME", "sgrnas"))
    orig_sgrna_list_file = config.get("orig_sgrna_list_file", config.get("ORIG_SGRNA_LIST_FILE",None))
    
    # Dual sgRNA support
    orig_sgrna2_list_file = config.get("orig_sgrna2_list_file", None)
    sgrna2_list_name = config.get("sgrna2_list_name", None)
    sgrna_pairs_file = config.get("sgrna_pairs_file", None)
    
    # Determine if this is a dual sgRNA library
    is_dual_sgrna = bool(orig_sgrna2_list_file and sgrna2_list_name)

    ignore_r2 = config.get("ignore_r2", "FALSE")
    trim5_lengths = config.get("trim5_lengths", "")
    keep_tmp = config.get("keep_tmp", "FALSE")
    save_unmapped = config.get("save_unmapped", "FALSE")
    
    # ReCo specific parameters
    vector_file = config.get("vector_file", "")
    expected_reads = config.get("expected_reads", 1000000)

    ncpu = config.get("ncpu", config.get("NCPU", 1))

    r1_seqs = config.get("r1_seqs")
    r2_seqs = config.get("r2_seqs")
    len_pattern = config.get("len_pattern", config.get("LEN_PATTERN", None))

    # Determine pattern for distance analysis (used by both reorientation and QC)
    pattern_to_use = None
    if config.get("len_pattern"):
        # Check if len_pattern is valid (only contains ACTG characters)
        if all(c.upper() in 'ACTG' for c in config['len_pattern']):
            pattern_to_use = config['len_pattern']

    # Convert string values to booleans if necessary
    if isinstance(search_revcomp, str):
        search_revcomp = search_revcomp.upper() == "TRUE"
    if isinstance(reformat_sgrna_list, str):
        reformat_sgrna_list = reformat_sgrna_list.upper() == "TRUE"
    if isinstance(skip_references, str):
        skip_references = skip_references.upper() == "TRUE"
    if isinstance(force_references, str):
        force_references = force_references.upper() == "TRUE"
    if isinstance(skip_fastq_processing, str):
        skip_fastq_processing = skip_fastq_processing.upper() == "TRUE"
    if isinstance(skip_matrix, str):
        skip_matrix = skip_matrix.upper() == "TRUE"
    if isinstance(skip_qc, str):
        skip_qc = skip_qc.upper() == "TRUE"
    if isinstance(skip_analysis, str):
        skip_analysis = skip_analysis.upper() == "TRUE"
    if isinstance(analysis_only, str):
        analysis_only = analysis_only.upper() == "TRUE"
    if isinstance(paired_end_fastq, str):
        paired_end_fastq = paired_end_fastq.upper() == "TRUE"
    if isinstance(keep_tmp, str):
        keep_tmp = keep_tmp.upper() == "TRUE"
    if isinstance(save_unmapped, str):
        save_unmapped = save_unmapped.upper() == "TRUE"

    generate_references_cmd_opts = ""
    
    config_file = config.get("config_file", "config.sh")

    global docker_paths

    if "docker_paths" in config or "DOCKER_PATHS" in config:
        docker_paths = config.get("docker_paths", config.get("DOCKER_PATHS"))
        if isinstance(docker_paths, str):
            docker_paths = docker_paths.upper() == "TRUE"
        logging.warning(f"DOCKER_PATHS is set to {docker_paths} in {config_file}, overriding any autodetection of docker.")
    else:
        docker_paths = util.is_docker_env()
    
    logging.info(f"docker_paths: {docker_paths}")

    if docker_paths:
        working_dir = config.get("docker_working_dir", config.get("DOCKER_WORKING_DIR"))
        fastq_dir = config.get("docker_fastq_dir", config.get("DOCKER_FASTQ_DIR"))
        output_dir = config.get("docker_output_dir", config.get("DOCKER_OUTPUT_DIR"))
    else:
        working_dir = config.get("working_dir", config.get("WORKING_DIR"))
        fastq_dir = config.get("fastq_dir", config.get("FASTQ_DIR"))
        output_dir = config.get("output_dir", config.get("OUTPUT_DIR"))

    # Auto-detect metadata file format (prefer YAML, fallback to TXT)
    # Do this after Docker paths are set up so we look in the right directory
    metadata_file = config.get("metadata_file", config.get("METADATA_FILE", None))
    if metadata_file is None:
        if os.path.exists("sample_metadata.yaml"):
            metadata_file = "sample_metadata.yaml"
        elif os.path.exists("sample_metadata.txt"):
            metadata_file = "sample_metadata.txt"
        else:
            logging.error("Neither sample_metadata.yaml nor sample_metadata.txt found.")
            logging.error("Please create a metadata file or specify one in the config.")
            sys.exit(1)
    
    logging.info(f"Using metadata file: {metadata_file}")

    # set default references dir
    references_dir = config.get("references_dir", config.get("REFERENCES_DIR", f"{output_dir}/references"))

    # Create references directory if it does not exist
    if not os.path.exists(references_dir):
        os.makedirs(references_dir)
        logging.info(f"Created references directory: {references_dir}")

    # process metadata
    libraries, sample_names, fastq_files = pm.process_metadata(metadata_file, fastq_dir)

    logging.info(f"working_dir: {working_dir}")
    logging.info(f"fastq_dir: {fastq_dir}")
    logging.info(f"output_dir: {output_dir}")
    logging.info(f"mode: {mode}")
    logging.info(f"count_methods: {count_methods}")
    logging.info(f"search revcomp: {search_revcomp}")
    logging.info(f"ncpu: {ncpu}")
    
    # Log dual sgRNA configuration
    if is_dual_sgrna:
        if "mageck" in count_methods:
            logging.info("Dual sgRNA library detected - MAGeCK2 will be used for counting")
        if "reco" in count_methods:
            logging.warning("Dual sgRNA libraries are not fully supported with ReCo. Using single library ReCo approach.")
        logging.info(f"Primary library: {sgrna_list_name} (from {orig_sgrna_list_file})")
        logging.info(f"Secondary library: {sgrna2_list_name} (from {orig_sgrna2_list_file})")
        
        # Validate tool availability based on count methods
        if "mageck" in count_methods:
            # Validate MAGeCK2 availability
            if not util.run_command("which mageck2"):
                logging.error("MAGeCK2 is required for dual sgRNA libraries but is not installed or not in PATH.")
                logging.error("Please install MAGeCK2 or use single sgRNA library configuration.")
                sys.exit(1)
        if "reco" in count_methods:
            # ReCo is only available for fastq mode
            if mode != "fastq":
                logging.error("ReCo only supports 'fastq' mode. Dual sgRNA libraries with ReCo require fastq mode.")
                sys.exit(1)
            # Check ReCo availability
            try:
                import reco
                logging.info("ReCo package found")
            except ImportError:
                logging.error("ReCo package is required but not installed. Please install ReCo.")
                sys.exit(1)
        
        # Ensure fastq mode for dual sgRNA
        if mode != "fastq":
            logging.error(f"Mode '{mode}' is not supported for dual sgRNA libraries. Only 'fastq' mode is supported.")
            sys.exit(1)
    else:
        if "mageck" in count_methods:
            logging.info("Single sgRNA library - MAGeCK will be used for counting")
        if "reco" in count_methods:
            # ReCo is only available for fastq mode
            if mode != "fastq":
                logging.error("ReCo only supports 'fastq' mode. Please change mode to 'fastq' or use 'mageck' count_method.")
                sys.exit(1)
            logging.info("Single sgRNA library - ReCo will be used for counting")
            # Check ReCo availability
            try:
                import reco
                logging.info("ReCo package found")
            except ImportError:
                logging.error("ReCo package is required but not installed. Please install ReCo.")
                sys.exit(1)

    parser = argparse.ArgumentParser(description="Run CRISPR screen pipeline")
    parser.add_argument("--skip-matrix", action="store_true", default=skip_matrix, help="Skip matrix generation")
    parser.add_argument("--skip-qc", action="store_true", default=skip_qc, help="Skip QC section")
    parser.add_argument("--skip-analysis", action="store_true", default=skip_analysis, help="Skip analysis section")
    parser.add_argument("--analysis-only", action="store_true", default=analysis_only, help="Run analysis section only")
    parser.add_argument("--skip-fastq-processing", action="store_true", default=skip_fastq_processing, help="Skip FASTQ processing")
    parser.add_argument("--skip-references", action="store_true", default=skip_references, help="Skip reference generation")
    parser.add_argument("--force-references", action="store_true", default=force_references, help="Force reference generation")
    parser.add_argument("--reformat-sgrna-list", action="store_true", default=reformat_sgrna_list, help="Reformat sgRNA list")
    parser.add_argument("-q", "--quiet", action="store_true", help="Suppress logging output")

    args = parser.parse_args()

    global quiet
    quiet = args.quiet

    skip_matrix = args.skip_matrix
    skip_qc = args.skip_qc
    skip_analysis = args.skip_analysis
    analysis_only = args.analysis_only
    skip_fastq_processing = args.skip_fastq_processing
    skip_references = args.skip_references
    force_references = args.force_references
    reformat_sgrna_list = args.reformat_sgrna_list

    generate_references_cmd_opts = ""
    if reformat_sgrna_list:
        generate_references_cmd_opts += " --reformat-sgrna-list"

    if analysis_only:
        logging.info("Running analysis section only.")
        logging.info(f"Config file: {config_file}")
        ma.mageck_test(config_file=config_file)
        sys.exit(0)


    if mode == "bam":
        required_software = ["bowtie2", "cutadapt", "samtools"]
        for software in required_software:
            if not util.run_command(f"which {software}"):
                logging.error(f"required software is not installed or on the current PATH: {software}")
                sys.exit(1)

        if skip_references:
            logging.info("skipping reference generation.")
        else:
            if os.path.exists(f"{references_dir}/{sgrna_list_name}.fa") and not force_references:
                logging.info("file exists, skipping reference generation.")
            else:
                output = util.run_command(f". {working_dir}/generate_reference_files.sh {generate_references_cmd_opts}")
                logging.info(output)

            if os.path.exists(f"{references_dir}/{sgrna_list_name}.1.bt2") and not force_references:
                logging.info("file exists, skipping bt2 index generation.")
            else:
                output = util.run_command(f". {working_dir}/make_bt2_index.sh")
                logging.info(output)
    # process each sample
    for i, library in enumerate(libraries):
        sample = sample_names[i]
        logging.info(f"{library}: {sample}")

        if not sample or config.get("force_library_as_name", "FALSE").upper() == "TRUE":
            if config.get("force_library_as_name", "FALSE").upper() != "TRUE" and config.get("force_library_as_name", "FALSE").upper() != "FALSE":
                logging.warning(f"Invalid value for force_library_as_name: {config.get('force_library_as_name')}. Assuming default False.")
            sample = library

        in_r1_fastq = fastq_files[i].split(";")[0]
        in_r2_fastq = fastq_files[i].split(";")[1]

        r1_space = in_r1_fastq.replace(",", " ")
        r2_space = in_r2_fastq.replace(",", " ")

        if not skip_fastq_processing:
            for file in r1_space.split() + r2_space.split():
                if not os.path.exists(file) or os.path.getsize(file) == 0:
                    logging.warning(f"{file} does not exist or is empty. Skipping library {library}.")
                    continue

        in_r1_fastq = f"{output_dir}/{library}_combined_R1.fq.gz"
        in_r2_fastq = f"{output_dir}/{library}_combined_R2.fq.gz"
        in_r1_re_fastq = f"{output_dir}/{library}_combined_R1.reoriented.fq.gz"
        in_r2_re_fastq = f"{output_dir}/{library}_combined_R2.reoriented.fq.gz"

        if skip_fastq_processing:
            logging.info("Skipping FASTQ processing.")
        else:
            if not os.path.exists(in_r1_fastq) or os.path.getsize(in_r1_fastq) == 0:
                logging.info(f"merging {r1_space}")
                logging.info(f"to: {in_r1_fastq}")
                util.run_command(f"cat {r1_space} > {in_r1_fastq}")
            else:
                logging.info(f"Merged file {in_r1_fastq} exists.")

            if paired_end_fastq:
                if not os.path.exists(in_r2_fastq) or os.path.getsize(in_r2_fastq) == 0:
                    logging.info(f"merging {r2_space}")
                    logging.info(f"to: {in_r2_fastq}")
                    util.run_command(f"cat {r2_space} > {in_r2_fastq}")
                else:
                    logging.info(f"Merged file {in_r2_fastq} exists.")
            else:
                logging.warning("skipping R2 because paired_end_fastq is set to false.")

            fastq_stats_file = f"{output_dir}/fastq_stats.txt"

            if not os.path.exists(fastq_stats_file):
                with open(fastq_stats_file, "w") as f:
                    f.write("library\tread_count\n")

            read_count = -1
            with open(fastq_stats_file) as f:
                for line in f:
                    if line.startswith(library):
                        read_count = int(line.split()[1])
                        break

            if read_count <= 0:
                logging.info(f"Computing read count for {library}")
                new_read_count = int(util.run_command(f"zcat {in_r1_fastq} | wc -l")) // 4
                logging.info(f"read count: {new_read_count}")

                with open(fastq_stats_file, "a") as f:
                    f.write(f"{library}\t{new_read_count}\n")
           
            # If no valid len_pattern, try to use first r1_seq
            if not pattern_to_use and r1_seqs:
                first_r1_seq = r1_seqs[0]
                if all(c.upper() in 'ACTG' for c in first_r1_seq):
                    pattern_to_use = first_r1_seq

            if reorient_fastq and mode != "bam":
                if config.get("paired_end_fastq") == "FALSE":
                    logging.warning("skipping reorientation because paired_end_fastq is set to false.")
                else:
                    if not os.path.exists(in_r1_re_fastq) or not os.path.exists(in_r2_re_fastq):
                        logging.info("reorienting fastq files (this will take a long time)...")
                        start_time = time.time()
                        # reorient_output = util.run_command(f"python {working_dir}/reorient_fastq_parallel.py --sequences-r1 {','.join(r1_seqs)} --sequences-r2 {','.join(r2_seqs)} --cpus {ncpu} --in-fastq-r1 {in_r1_fastq} --in-fastq-r2 {in_r2_fastq} --out-fastq-r1 {in_r1_re_fastq} --out-fastq-r2 {in_r2_re_fastq} --plot --plot-prefix {output_dir}/{sample}")
                        reorient_output = reorient.main(
                            sequences_r1=','.join(r1_seqs),
                            sequences_r2=','.join(r2_seqs),
                            cpus=ncpu,
                            in_fastq_r1=in_r1_fastq,
                            in_fastq_r2=in_r2_fastq,
                            out_fastq_r1=in_r1_re_fastq,
                            out_fastq_r2=in_r2_re_fastq,
                            plot=True,
                            plot_prefix=f"{output_dir}/{sample}",
                            plot_only=False,
                            validate=False
                        )
                        end_time = time.time()

                        logging.info(f"Reorientation completed in {end_time - start_time} seconds")
                        logging.info(reorient_output)
                    else:
                        logging.info("reoriented fastq files already exist, skipping reorientation.")

        if skip_qc:
            logging.warning("skipping QC section.")
        else:
            logging.info("--------")
            logging.info("   R1")
            logging.info("--------")
            # r1_top_seqs_output = util.run_command(f"python {working_dir}/find_top_sequences.py --in-fastq {in_r1_fastq} --out-fasta {output_dir}/{library}_R1.fa --out-alignment-clustalw {output_dir}/{library}_R1.clustalw.aln --sample-size 1000")
            r1_top_seqs_output = top.find_top_sequences(
                in_fastq=in_r1_fastq,
                out_fasta=f"{output_dir}/{library}_R1.fa",
                out_alignment_clustalw=f"{output_dir}/{library}_R1.clustalw.aln"
            )
            if r1_top_seqs_output:
                for line in r1_top_seqs_output.splitlines():
                    sys.stdout.write(line + '\n')

            # Run distance analysis if we have a valid pattern
            if pattern_to_use:
                r1_distance_to_pattern_output = dist.distance_to_pattern_frequencies(
                    in_fastq=in_r1_fastq,
                    pattern=pattern_to_use,
                    nreads=10000,
                    start_at=0,
                    stop_at=0,
                    center_sequence=True,
                    frac_threshold=0.005
                )
                if r1_distance_to_pattern_output:
                    for line in r1_distance_to_pattern_output.splitlines():
                        logging.info(line)

            if config.get("paired_end_fastq") == "FALSE":
                logging.warning("skipping R2 because paired_end_fastq is set to FALSE.")
            else:
                logging.info("--------")
                logging.info("   R2")
                logging.info("--------")
                r2_top_seqs_output = top.find_top_sequences(
                    in_fastq=in_r2_fastq,
                    out_fasta=f"{output_dir}/{library}_R2.fa",
                    out_alignment_clustalw=f"{output_dir}/{library}_R2.clustalw.aln"
                )
                if r2_top_seqs_output:
                    for line in r2_top_seqs_output.splitlines():
                        sys.stdout.write(line + '\n')

                if pattern_to_use:
                    r2_distance_to_pattern_output = dist.distance_to_pattern_frequencies(
                        in_fastq=in_r2_fastq,
                        pattern=pattern_to_use,
                        nreads=10000,
                        start_at=0,
                        stop_at=0,
                        center_sequence=True,
                        frac_threshold=0.005
                    )
                    if r2_distance_to_pattern_output:
                        for line in r2_distance_to_pattern_output.splitlines():
                            logging.info(line)

            if os.path.exists(in_r1_re_fastq):
                logging.info("-------------------")
                logging.info("   R1 reoriented")
                logging.info("-------------------")
                r1_reoriented_top_seqs_output = top.find_top_sequences(
                    in_fastq=in_r1_re_fastq,
                    out_fasta=f"{output_dir}/{library}_R1_reoriented.fa",
                    out_alignment_clustalw=f"{output_dir}/{library}_R1_reoriented.clustalw.aln"
                )
                if r1_reoriented_top_seqs_output:
                    for line in r1_reoriented_top_seqs_output.splitlines():
                        sys.stdout.write(line + '\n')

                if pattern_to_use:
                    r1_reoriented_distance_to_pattern_output = dist.distance_to_pattern_frequencies(
                        in_fastq=in_r1_re_fastq,
                        pattern=pattern_to_use,
                        nreads=10000,
                        start_at=0,
                        stop_at=0,
                        center_sequence=True,
                        frac_threshold=0.005
                    )
                    if r1_reoriented_distance_to_pattern_output:
                        for line in r1_reoriented_distance_to_pattern_output.splitlines():
                            logging.info(line)
            else:
                logging.warning("Reoriented R1 fastq files do not exist. Skipping reoriented QC.")

            if config.get("paired_end_fastq") == "FALSE":
                logging.warning("skipping R2 because paired_end_fastq is set to FALSE.")
            else:
                if os.path.exists(in_r2_re_fastq):
                    logging.info("-------------------")
                    logging.info("   R2 reoriented")
                    logging.info("-------------------")
                    r2_reoriented_top_seqs_output = top.find_top_sequences(
                        in_fastq=in_r2_re_fastq,
                        out_fasta=f"{output_dir}/{library}_R2_reoriented.fa",
                        out_alignment_clustalw=f"{output_dir}/{library}_R2_reoriented.clustalw.aln"
                    )
                    if r2_reoriented_top_seqs_output:
                        for line in r2_reoriented_top_seqs_output.splitlines():
                            sys.stdout.write(line + '\n')

                    if pattern_to_use:
                        r2_reoriented_distance_to_pattern_output = dist.distance_to_pattern_frequencies(
                            in_fastq=in_r2_re_fastq,
                            pattern=pattern_to_use,
                            nreads=10000,
                            start_at=0,
                            stop_at=0,
                            center_sequence=True,
                            frac_threshold=0.005
                        )
                        if r2_reoriented_distance_to_pattern_output:
                            for line in r2_reoriented_distance_to_pattern_output.splitlines():
                                logging.info(line)
                else:
                    logging.warning("Reoriented R2 fastq files do not exist. Skipping reoriented QC.")

        if reorient_fastq and mode != "bam" and paired_end_fastq:
            in_r1_fastq = in_r1_re_fastq
            in_r2_fastq = in_r2_re_fastq

        logging.info(f"r1: {in_r1_fastq}")
        logging.info(f"r2: {in_r2_fastq}")

        if mode == "grep":
            count_output = cg.count_grep(
                sample=sample,
                mode=mode,
                output_dir=output_dir,
                in_r1_fastq=in_r1_fastq,
                in_r2_fastq=in_r2_fastq,
                sgrna_list_name=sgrna_list_name,
                ncpu=ncpu
            )
        elif mode == "bam":
            # BAM mode is not supported for dual sgRNA libraries
            if is_dual_sgrna:
                logging.error("BAM mode is not supported for dual sgRNA libraries. Please use 'fastq' mode.")
                sys.exit(1)
            
            # BAM mode always uses MAGeCK (ReCo only supports fastq mode)
            if "reco" in count_methods:
                logging.info("BAM mode always uses MAGeCK for counting. ReCo is only available for fastq mode.")
                
            required_software = ["bowtie2", "cutadapt", "samtools"]
            for software in required_software:
                if not util.run_command(f"which {software}"):
                    logging.error(f"required software is not installed or on the current PATH: {software}")
                    sys.exit(1)
            trim_output = util.run_command(f". {working_dir}/trim_reads.sh")
            logging.info(trim_output)
            align_output = util.run_command(f". {working_dir}/align_bt2.sh")
            logging.info(align_output)

            # Always use MAGeCK for BAM mode
            count_output = count.run_mageck_count(
                mode=mode,
                sample=sample,
                output_dir=output_dir,
                references_dir=references_dir,
                sgrna_list_name=sgrna_list_name,
                in_r1_fastq=in_r1_fastq,
                in_r2_fastq=in_r2_fastq,
                trim5_lengths=trim5_lengths,
                keep_tmp=keep_tmp,
                save_unmapped=save_unmapped,
                search_revcomp=search_revcomp,
                ignore_r2=ignore_r2
            )           
            logging.info(count_output)
        elif mode == "fastq":
            # Run counting for each specified method
            for count_method in count_methods:
                logging.info(f"Running counting with method: {count_method}")
                
                if is_dual_sgrna:
                    if count_method == "mageck":
                        # Use MAGeCK2 for dual sgRNA libraries
                        count_output = count2.run_mageck2_count(
                            mode=mode,
                            sample=sample,
                            output_dir=output_dir,
                            references_dir=references_dir,
                            sgrna_list_name=sgrna_list_name,
                            sgrna2_list_name=sgrna2_list_name,
                            in_r1_fastq=in_r1_fastq,
                            in_r2_fastq=in_r2_fastq,
                            sgrna_pairs_file=sgrna_pairs_file,
                            trim5_lengths=trim5_lengths,
                            keep_tmp=keep_tmp,
                            save_unmapped=save_unmapped,
                            search_revcomp=search_revcomp,
                            search_revcomp2=search_revcomp,  # For now, use same setting for both libraries
                            ignore_r2=ignore_r2
                        )
                    elif count_method == "reco":
                        # Use ReCo for dual sgRNA libraries (single library approach)
                        logging.warning("Using single library ReCo approach for dual sgRNA library")
                        count_output = count_reco.run_reco_count(
                            mode=mode,
                            sample=sample,
                            output_dir=output_dir,
                            references_dir=references_dir,
                            sgrna_list_name=sgrna_list_name,
                            in_r1_fastq=in_r1_fastq,
                            in_r2_fastq=in_r2_fastq,
                            trim5_lengths=trim5_lengths,
                            keep_tmp=keep_tmp,
                            save_unmapped=save_unmapped,
                            search_revcomp=search_revcomp,
                            ignore_r2=ignore_r2,
                            cores=ncpu,
                            vector_file=vector_file,
                            expected_reads=expected_reads
                        )
                else:
                    if count_method == "mageck":
                        # Use regular MAGeCK for single sgRNA libraries
                        count_output = count.run_mageck_count(
                            mode=mode,
                            sample=sample,
                            output_dir=output_dir,
                            references_dir=references_dir,
                            sgrna_list_name=sgrna_list_name,
                            in_r1_fastq=in_r1_fastq,
                            in_r2_fastq=in_r2_fastq,
                            trim5_lengths=trim5_lengths,
                            keep_tmp=keep_tmp,
                            save_unmapped=save_unmapped,
                            search_revcomp=search_revcomp,
                            ignore_r2=ignore_r2
                        )
                    elif count_method == "reco":
                        # Use ReCo for single sgRNA libraries
                        count_output = count_reco.run_reco_count(
                            mode=mode,
                            sample=sample,
                            output_dir=output_dir,
                            references_dir=references_dir,
                            sgrna_list_name=sgrna_list_name,
                            in_r1_fastq=in_r1_fastq,
                            in_r2_fastq=in_r2_fastq,
                            trim5_lengths=trim5_lengths,
                            keep_tmp=keep_tmp,
                            save_unmapped=save_unmapped,
                            search_revcomp=search_revcomp,
                            ignore_r2=ignore_r2,
                            cores=ncpu,
                            vector_file=vector_file,
                            expected_reads=expected_reads
                        )
                
                logging.info(f"Completed counting with {count_method}: {count_output}")
        else:
            logging.error(f"Invalid mode: {mode}")
            sys.exit(1)

    if skip_matrix:
        logging.info("Skipping matrix generation.")
    else:
        logging.info("building count and cpm matrix files...")
        # For cpm_matrix.py - explicitly pass the detected metadata file
        cpm_matrix_output = cpm.count_and_cpm_matrix(config_file=config_file, metadata_file=metadata_file)
        if cpm_matrix_output:
            for line in cpm_matrix_output.splitlines():
                logging.info(line)
        logging.info("plotting count distribution histograms...")
        # For count_distribution_histograms.R - pass metadata file as environment variable
        histogram_output = util.run_command(f"METADATA_FILE={metadata_file} Rscript {working_dir}/count_distribution_histograms.R")
        if histogram_output:
            for line in histogram_output.splitlines():
                logging.info(line)

        if skip_analysis:
            logging.info("Skipping analysis section.")
        else:
            # Run analysis based on which count methods were used
            if "mageck" in count_methods:
                logging.info(f"Running MAGeCK analysis...")
                analysis_output = util.run_command(f". {working_dir}/mageck_analysis.sh --mode {mode}")
                if analysis_output:
                    for line in analysis_output.splitlines():
                        logging.info(line)
            
            if "reco" in count_methods:
                logging.info("ReCo analysis completed during counting step. No additional analysis needed.")
            
            # Check for unknown methods (shouldn't happen due to validation, but good to be safe)
            unknown_methods = [method for method in count_methods if method not in ["mageck", "reco"]]
            if unknown_methods:
                logging.warning(f"Unknown analysis methods: {unknown_methods}")
if __name__ == "__main__":
    main()
