# CRISPR screen sequence analysis
This analysis pipeline takes in sequence data (paired-end fastq files) containing sgRNA sequences, along with a list of known sgRNAs (the CRISPR library) and performs counting and downstream analysis using MaGeCK.

## purpose
* wrap configuration, settings, and tools around MaGeCK for counting and analysis of CRISPR screen sequencing data.
* inspect and correct input reads to be compatible with MaGeCK

A main feature of this repository is `reorient_fastq_parallel.py` which corrects fastq files where DNA fragments were sequenced in both forward and reverse directions, i.e. the "forward" end may be in either R1 or R2.   In its current form, MaGeCK assumes all (R1) reads originate from the same end of the amplicon.  While MaGeCK is capable of using R2 to aid in finding sgRNA sequences, it will search R2 only as the reverse complement of what it is searching for in R1.  A QC tool, `find_top_sequences.py`, is included to help determine whether read reorientation is necessary.

# software requirements (this list is not exhaustive; check out docker/env.yaml)
* mageck 0.5.9.5
* clustalw 2.1
* seqtk 1.4
* Python 3.9
* Python packages:
  * pyyaml
  * Biopython
  * termcolor
  * scikit-learn
  * matplotlib
  * scipy
  * numpy
* R 3.4
* R packages
  * tidyverse
  * openxlsx
  * org.Hs.eg.db
  * ggrepel
  * GGally
* system/cmd line utilities
  * column (from util-linux package)
  * gawk (awk pretending to be gawk will work)
  * file
  * jq
  * yq
  * cmake
* alignment/BAM mode requirements (optional for fastq mode)
  * bowtie2 2.4.5 (2.5.2 ok)
  * samtools 1.6 (1.18 ok)
  * cutadapt 4.1 (4.5 ok)

## optional software
* pigz
* GNU parallel (parsort)

## conda setup
Conda is not required, but setting up a conda environment may be convenient.  You could use `docker/env.yaml`
1. install miniconda, e.g. [https://docs.conda.io/en/latest/miniconda.html#linux-installers](https://docs.conda.io/en/latest/miniconda.html#linux-installers)
2. from a fresh miniconda install:
  ```
  conda init
  conda create --name crispr --f env.yaml
  ```

## docker
A docker option is outlined below.  This repo includes two files, `docker/Dockerfile` and `docker/env.yaml`.  The micromamba approach is just one way to build an image.  Experienced users may prefer to do it a different way.

### 1. pull the docker image and note the image name
The below command will pull version 'latest'
  ```
  docker pull mambaorg/micromamba
  ```
### 2. build a new docker image
It will which will contain all the software needed as well as ensure the environment is activated when a docker container is run. tag as desired
  ```
  cd docker
  docker build --quiet --tag my_pipelines/crispr:0.4 .
  cd -
  ```
### 4. run the docker container
  ```
  $ base_dir=$(pwd)
  $ docker run --rm -idt \
   --name crispr-pipeline \
   -v $base_dir:/repo \
   -v $base_dir/output:/output \
   -v $base_dir/input:/input \
  my_pipelines/crispr:0.1 \
  bash
  $ docker attach crispr-pipeline
  (base) mambauser@397fd3feaf7b:/repo$ python run_pipeline.py
  ```

# pipeline setup

### 1. acquire sgRNA table for the given CRISPR library.
* Place it in the working directory (repo base directory). Name it`sgRNAs.txt`, for example.

* The file should be tab-delimited plain text with Unix style line endings (LF and *not* CRLF).

* The file should be formatted for MAGeCK. That is, it should have 3 columns with a header: sgRNA_ID, sgRNA_sequence, target_name.  The header text can be anything.

* Optional: use the validation tool: `bash validate_sgrna_list.sh --file sgRNAs.txt`

### 2. create the sample metadata file.
* Place in base repo dir and name it `sample_metadata.txt`.
* Use the example as a template `examples/sample_metadata.txt`.
* The metadata file should be a 4 column tab-delimited text file with header.  The column names do not matter, but the order does:
  * **column 1** - library
A Library ID.
  * **column 2** - sample
A sample name/alias.  Short and readable withno spaces or special characters.  
  * **column 3** - fastq_r1
A comma-separated list of fastq files for R1 (e.g. `library_lane1_R1.fq.gz,library_lane2_R1.fq.gz`). Full path or file name only are accepted. If file name only, fastq path (e.g. `input/`) must be set in `config.yaml`.  
  * **column 4** - fastq_r2
A comma-separated list of fastq files for R2 (e.g. `library_lane1_R2.fq.gz,library_lane2_R2.fq.gz`).

### 3. input/output dirs
* make the dirs: 
  * `mkdir input`
  * `mkdir output`
* copy fastq files into `input/`

### 4. create configuration script called `config.yaml`.  Edit as needed. 
* Use the example as a template `examples/config.yaml`. 
* Note: if running within a docker container, be sure to set docker-specific paths and options.
* ensure input/output dirs are specified.

### 5. run pipeline

* run the pipeline using: `python run_pipeline.py` in the repo directory. 
* refer to MAGeCK documentation for interpreting the output of `mageck count` or `mageck test`.  


# run specific tools

## 1. Find top sequences with `find_top_sequences.py`

### A. The minimum input parameter required is a fastq file, e.g.

```
python find_top_sequences.py --in-fastq sample.fastq.gz
```
  
* The output will consist of the top sequence classes, pruned at threshold 0.1 (must represent > 10% of the reads).
* The output will also print a degenerate consensus sequence derived from a multiple sequence alignment (ClustalW)
* It is recommended to run both R1 and R2 fastq files.
* examine all parameters with `python find_top_sequences.py --help`

### B. Sampling

* Sampling is performed by seqtk. By default, seqtk uses the same random seed, so results are reproducible. This seed can be explicitly set using `--seed`.

### C. Classification

* The number of reads sampled for classification can be set using `--sample-size` (Default: 1000).
* The default classification method uses SequenceMatcher from difflib.
* A second method can be used by setting `--method-2`.  A model must also be specified using `--model <model>`.  An included model is `crispr_fq_model4c`.
* Method 2 is a multinomial naive Bayes classifier (using MultinomialNB from scikit-learn). It can be much faster than the default method for large n. The downside to this method is the training set must use reads from the same sgRNA vector as the test set.  Essentially, prior knowledge of and access to sequences for training is required. A training set (`read_training_set.txt`) and model (`crispr_fq_model4c`) are provided, but there is no guarantee they work with every sytem. A new model can be prepared using `find_top_sequences_train.py`.
* Note: when running Method 2, it reports OOD (Out Of Distribution) sequences. These sequences are essentially "unseen" or "novel" data to the model.  These could be used to retrain and improve the model, especially if they represent a high percentage of samples reads.  

### D. Multiple sequence alignment

* MSA related output files are saved by default to the base dir of the input fastq file. They can be specified with `--out-fasta` and `--out-alignment-clustalw`. 
* MSA is very slow for large sequence sets. Therefore, by default, the MSA sample size ceiling is set to 200, even if `--sample-size` is set much higher. This ceiling can be set using `--msa-sample-size`. Use caution when increasing this value, as the MSA can become very slow.
* The MSA is mainly used to generate a consensus and as such gaps are not desired, so they are heavily penalized. 
* The resulting consensus sequence is degenerate and uses [IUPAC notation](https://en.wikipedia.org/wiki/Nucleic_acid_notation).

### E. Interpretation

* Read pairs correctly oriented:
  * There will be one dominant sequence class.
  * The consensus sequence will be well-defined throughout, the notable exception being the position of the sgRNA, which is highly variable. Here will be 20 bases of Ns.
Number of sequence classes: 3
```
Class 1:
Representative sequence: TTGTGGAAAGGACGAAACACCGTAGAAGAGGATGTCAAACGTGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTGAATTCGCTAGCTAGGTCTTGAAAGGA
Number of sequences in class: 569 (56.90%)
Class 2:
Representative sequence: CCAATTCCCACTCCTTTCAAGACCTAGCTAGCGAATTCAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAACCCTTCCTGTATTGTCTCATTCGGTGTTTCGT
Number of sequences in class: 430 (43.00%)
Class 3:
Representative sequence: CCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCAAAAATTAAAAAAAATAAAAAAAAAAAAAAAAAATAAAAAAAAAAATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTAT
Number of sequences in class: 1 (0.10%)
```

![consensus for correctly oriented reads](https://github.com/mctp/crispr-pipeline/blob/main/consensus_correct.png?raw=true)

* Read pairs incorrectly oriented:
  * There will be 2 dominant sequence classes, one where R1 is the "forward" end of the DNA and one where R1 is the "reverse" end of the DNA.
  * The MSA consensus uses degenerate nucleotides (IUPAC) and will show poor consensus if the fastq is randomly oriented.  Stretches of Ns at both ends are likely due to sgRNA sequences.  

```
Number of OOD sequences: 1 (0.10%)
OOD sequences (first 20):
ATTCGTGGGTGTGTTTGGTAGCTTGTGGGTGTGTTTGGTAGCTGGGGGCCCTTGTTGGTAGCTTGTGGGTGTATTTGGTAGCTTGTGGGTGTGTTTTGAAGCTTGTGGGTGTGTTTGGTAGCTTGTGGGTGTGTTTGGTAGCTTGTGTGGG
Number of sequence classes: 1
Class 1:
Representative sequence: CCAATTCCCACTCCTTTCAAGACCTAGCTAGCGAATTCAAAAAAGCACCGACTCGGTGCCACTTTTTCAAGTTGATAACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAACCAGTAAGCTTAGGTGGGCAGCGGTGTTTCGT
Number of sequences in class: 999 (99.90%)
```

![consensus for incorrectly oriented reads](https://github.com/mctp/crispr-pipeline/blob/main/consensus_incorrect.png?raw=true)


## 2. Reorient fastq files with `reorient_fastq_parallel.py`
### A. Minimum required parameters.

5 fastq files: The input R1 and R2, the output R1 and R2, and the output interleaved file

At least 1 search sequence (--sequences-r1 must not be empty).  

```
python reorient_fastq_files.py \ 
  --in-fastq-r1 sample_R1.fq.gz \
  --in-fastq-r2 sample_R2.fq.gz \
  --out-fastq-r1 sample_R1.reoriented.fq.gz \
  --out-fastq-r2 sample_R2.reoriented.fq.gz \
  --file-interleaved-gz interleaved.fq.gz \
  --sequences-r1 <sequence1>,<sequence2> \
  --sequences-r2 <sequence3>,<sequence4>
```
The search sequence(s) must be known in advance.  Choose a sequence expected to be found in R1 for properly oriented read pairs.  This sequence should be long enough to be unique, say 18-20 bases, but not so long that exact matching is thwarted by sequencing errors.
The example `config.sh` file contains this snippet with 2 search sequences for R1 and R2:
```
# expected sequences in R1
r1_seqs=(GCTATTTCTAGCTCTAAAAC                 TCCCACTCCTTTCAAGA)
#        ^^^ left adjacent to sgRNA           ^^^ PCR primer (rev) part (5' 5 bases removed)

# expected sequences in R2
r2_seqs=(GTTTTAGAGCTAGAAATAGC                 GAAAGGACGAAACACCG)
#        ^^^ left adjacent to sgRNA (revcomp) ^^^ PCR primer (fwd) part (5' 5 bases removed)
```
### B. Optional parameters.
* Use `-p|--cpus` to take advantage of multiple processors. 
* Batch size for processing can be changed from default (1000) using `--batch-size` but changing this value is not recommended.
* By default, all reads are processed.  To specify a number of reads, use `--nreads`. 
* `--validate` is used to confirm that the read IDs in the input R1 and R2 fastq files are found in the output R1 and R2 files.    

### C. Plots
* If you want to visualize the effect of read reorientation (before and after), use `--plot`.  To only do the plot (reorientation was previously performed) also add `--plot-only`.
* Plotting is performed by `read_reorientation_barplot.py`.  It outputs a PNG and PDF of the actual plot as well as a table of results.  To specify a prefix for these output files, use `--plot-prefix`. 
* Briefly, reads are sampled and the sequences specified by `--sequences-r1` are searched for.  Hits are tallied and a stacked barplot for each fastq file is produced showing the percent of reads sampled containing the query sequences (Found) or not (Not Found) is shown.  The sequence choice can be overridden by using `--plot-search-sequence`.

# downstream analysis / testing
### building counts and CPM matrices
* The script `cpm_matrix.py` runs by default with `run_pipeline.py` but it can be run standalone, also, either specifying the configuration file
```
python cpm_matrix.py --config config.sh
```
or by specifying the sample metadata file and input fastq dir:
```
python cpm_matrix.py --metadata sample_metadata.txt --input-fastq-dir input
```
Counts are obtained for all samples in the metadata file not commented out.

### plot count histograms
* After running `cpm_matrix.py`, histograms are produced by `run_pipeline.sh` using `count_distribution_histograms.R`.
* To run separately, the pipeline context is still required (config,sh, sample_metadata.txt).
```
Rscript count_distribution_histogram.R
```
* for each sample, 4 plots are produced for counts, log10 counts, CPM, and log10 CPM.  Plots containing all samples are also produced (counts, log10 coutns, CPM, and log10 CPM).
* The terminal output also prints count stats (from mageck count summary files) and this also includes coverage of mapped reads.

### run/skip the analysis
* by default, `mageck test` is run using comparison details in `config.sh`.  If a comparison is not configured, the pipeline will run correctly and generate counts, but may exit with an error.
* use `--skip-analysis` to stop processing data after running `mageck count`.  
* use `--analysis-only` to perform the analysis without processing fastq and counting.  Often the counting process only needs to be done one time, but comparisons may be tried in different ways, depending on the experiment.
* Comparisons are done with `mageck_analysis.sh` and currently only runs `mageck test` with options specified in the configuration file.  

### misc
* use `--help` with any python script to view command line parameters.
* `run_pipeline.py` and `config.yaml` are superseding `run_pipeline.sh` and `config.sh` (respectively).
* optimal settings for `reorient_fastq_parallel.py` are `--batch-size 1000` and `--chunk-size 1000000`.  `--cpus` should be set to something >1 where possible, such as 4, 8 or 16; using higher does not help much.  1M chunk size keeps memory overhead low.  these numbers were obtained empirically and so YMMV.  
