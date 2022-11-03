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
* count reads


# software requirements
bowtie2 2.4.5
samtools 1.6
cutadapt 4.1
mageck 0.5.9.5
a linux-like command-line environment

## optional software
GNU parallel (we use parsort if available)

## conda setup
conda is not required, but setting up a conda environment may be convenient.
1. install miniconda, e.g. [https://docs.conda.io/en/latest/miniconda.html#linux-installers](https://docs.conda.io/en/latest/miniconda.html#linux-installers)
2. from a fresh miniconda install:
  ```
  conda init
  conda create -name crispr
  conda install -c bioconda bowtie2=2.4.5 samtools=1.6 cutadapt=4.1 mageck=0.5.9.5
  conda install parallel
  ```

# pipeline setup

1. acquire sgRNA list file for CRISPR library. This code was written around a 2 column file (shRNA ID, sequence), which is reformatted by this code to work with MAGeCK.
2. create a sample metadata file (e.g. `sample_metadata.txt`).  This should be a 3 column tab-delimited text file with header. Column names do not matter, but order and format does.  
  a. col 1 (library) is a library ID  
  b. col 2 (sample) is a sample name/alias, preferably something easy to read.  But use no spaces and keep it short!  
  c. col 3 (fastqs) is a comma-separated list of fastq files (e.g. `R1.fq.fz,R2.fq.gz`). Full path or file name accepted (as long as the path to fastqs is specified in `config.sh`).  Currently, only 2 fastq files are allowed, R1 and R2, so combine in advance where R1 and/or R2 have multiple associated fastq files. Todo: remove this constraint.  
3. edit the `config.sh` file. Use the example values as a guide and heed the comments.

# run pipeline

1. run the pipeline using: `bash run_pipeline.sh` in the repo directory.
2. reference files and bowtie2 index need only be created once per CRISPR library.  After successful first run, manually comment out these lines from run_pipeline.sh:
  ```
  # source generate_reference_files.sh
  # source make_bt2_index.sh
  ```  
3. refer to MAGeCK documentation for interpreting the output.  
