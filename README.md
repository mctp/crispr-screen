# CRISPR screen sequence analysis
this small pipeline takes in sequence data (paired-end fastq files) containing sgRNA sequences, along with a list of known sgRNAs (the CRISPR library), and processes, aligns, and counts the reads.

## prevent under-counting caused by the following
* DNA fragments sequenced in forward and reverse direction
* variable position of sgRNA sequences within reads
* sequence mismatches due to PCR or sequencing error

## overview of steps
* generate references
  * format sgRNA library file
  * generate  a fasta file from sgRNA library containing flanking vector sequences
  * build aligner index
* trim reads
  * trim R1 and R2 separately using 5' vector sequence adjacent to sgRNA
  * return only matching (trimmed) reads, cut to 20 bases and not less than 18 bases
  * combine R1 and R2 results into single fastq file
* align modified reads to custom reference index
* count reads

# software requirements
* bowtie2 2.4.5 (2.5.2 ok)
* samtools 1.6 (1.18 ok)
* cutadapt 4.1 (4.5 ok)
* mageck 0.5.9.5
* a linux-like command-line environment with the following command utilities
  * column (from util-linux package)
  * gawk (awk pretending to be gawk will work)
  * file
  * dos2unix

## optional software
* GNU parallel (we use parsort if available)

## conda setup
conda is not required, but setting up a conda environment may be convenient.
1. install miniconda, e.g. [https://docs.conda.io/en/latest/miniconda.html#linux-installers](https://docs.conda.io/en/latest/miniconda.html#linux-installers)
2. from a fresh miniconda install:
  ```
  conda init
  conda create --name crispr -c bioconda -c conda-forge
  # required
  conda install bowtie2=2.4.5 samtools=1.6 cutadapt=4.1 mageck=0.5.9.5
  # optional
  conda install pigz parallel
  ```

## docker container suggestion
a docker option is outlined below.  this repo includes two files, `docker/Dockerfile` and `docker/env.yaml`, as suggestions.  docker option is not fully tested.
1. pull the docker image and note the image name (the below command will pull version 'latest')
  ```
  docker pull mambaorg/micromamba
  ```
2. build a new docker image, which will contain all the software needed as well as ensure the environment is activated when a docker container is run. tag as desired
  ```
  cd docker
  docker build --quiet --tag my_pipelines/crispr:0.1 .
  cd -
  ```
4. run the docker container / start the pipeline with a single command
  ```
  base_dir=$(pwd)
  # container runs in background and removes itself when finished.  use -ti instead of -dt if you want to see the pipeline output in the current terminal
  docker run --rm -dt --name crispr-pipeline \
    -v $base_dir:/repo \
    -v $base_dir/output:/output \
    -v $base_dir/input:/input \
    my_pipelines/crispr:0.1 \
    bash -c "cd /repo && bash run_pipeline.sh"
  ```
3.  run the docker container interactively, providing access to the repo code, input directory, and output directory
  ```
  base_dir=($pwd)
  # make dirs if they do not exist
  mkdir -p "$base_dir/output"
  mkdir -p "$base_dir/input"
  # cp fastq files into input dir
  docker run -dt --name crispr-pipeline \
    -v $base_dir:/repo \
    -v $base_dir/output:/output \
    -v $base_dir/input:/input \
    my_pipelines/crispr:0.1
  # get a cmd line prompt inside the container, ready to use
  docker exec -ti crispr-pipeline bash
  # test
  samtools --version
  # run pipeline
  cd /repo
  bash run_pipeline.sh
  # do work
  # exit container e.g. cntrl+d and stop the container if done
  docker stop crispr-pipeline
  docker rm crispr-pipeline
  ```
# pipeline setup

1. acquire sgRNA list file for CRISPR library and place it in the base repo. example: `examples/human.sgRNAs.txt`
    1. your sgRNA list should be formatted for MAGeCK (sgRNA ID, sgRNA sequence, target/gene name)  
    2. uncommon/legacy: if your sgRNA library file contains 2 columns (id/sequence) where the id column is in the form sg<gene_name>_<numeric_id> use the `--reformat-sgrna-list` cmd line option when running `run_pipeline.sh` 

2. create sample metadata file and place in base repo dir. example: `examples/sample_metadata.txt`.  This should be a 3 column tab-delimited text file with header. Column names do not matter, but order and format does.
    1. col 1 (library) is a library ID  
    2. col 2 (sample) is a sample name/alias, preferably something easy to read.  But use no spaces and keep it short!  
    3. col 3 (fastqs) is a comma-separated list of fastq files (e.g. `R1.fq.fz,R2.fq.gz`). Full path or file name accepted (as long as the path to fastqs is specified in `config.sh`).  Currently, only 2 fastq files are allowed, R1 and R2, so combine in advance where R1 and/or R2 have multiple associated fastq files. Todo: remove this constraint.  

3. create configuration script called `config.sh`.  Edit as needed. 
    1. Example: `examples/config.sh`. 
    2. Example for docker: `examples/config_docker.sh`.

# run pipeline

1. run the pipeline using: `bash run_pipeline.sh` in the repo directory.
2. reference files and bowtie2 index need only be created once per CRISPR library.  After a successful first run, `run_pipeline.sh` will skip these steps if files exist. These steps can also be explicitly skipped using the `--skip-references` command line option.  Conversely, using `--force` will overwrite the reference/index files.
3. refer to MAGeCK documentation for interpreting the output.  
