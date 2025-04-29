#------------------------------------------------#
##################################################
####   CRISPR screen pipeline configuration   ####
##################################################
#------------------------------------------------#

MODE=fastq 
# required
# possible values:
#   - fastq (recommended)
#   - bam (perform alignment and count from bam file)
#   - grep (experimental use only)

SEARCH_REVCOMP=TRUE # fastq mode: when counting reads, search for reverse complement of sgrna sequence
					# bam mode: when generating aligner reference, use reverse complement of the sgrna sequence
                    # try changing this if mapping rate is low

# fastq mode only
reorient_fastq=TRUE # if PCR product is blunt end ligated and thus orientation of reads are randomly flipped
ignore_r2=TRUE      # set FALSE if not using paired-end reads OR if R2 cannot cover the sgRNA sequence
trim5_lengths= # comma-separated numbers or leave empty for default (recommended)
keep_tmp=FALSE  # FALSE recommended
save_unmapped=FALSE # FALSE recommended

# notes:
#       save_unmapped TRUE is incompatible with trim5_lengths with more than 1 value
#       sgrna length forced to 20 if trim5_lengths, keep_tmp, and/or save_unmapped are set

#######################
####   resources   ####
#######################

NCPU=32

###################
####   paths   ####
###################

DOCKER_PATHS=TRUE # TRUE if paths are relative to docker container, FALSE if paths are relative to host
skip_qc=FALSE # TRUE if QC has already been run and you want to skip it

# path to code base (this repo)
WORKING_DIR=.
DOCKER_WORKING_DIR=/repo

# put fastq files in FASTQ_DIR
# if metadata provides full path to fastq dir, this can be left blank
FASTQ_DIR=input
DOCKER_FASTQ_DIR=/input

OUTPUT_DIR=output # this dir must exist. if using docker, write permissions may need to be explicitly set
DOCKER_OUTPUT_DIR=/output # this dir must exist. if using docker, write permissions may need to be explicitly set

METADATA_FILE=sample_metadata.txt # see sample_metadata.txt and sample_metadata.yaml in examples/
                                               # for valid metadata files


#########################
####   cloud paths   ####
#########################

cloud_storage_command="gcloud storage" # e.g. 'gcloud storage' or 'gsutil'
cloud_storage_ls_command="gcloud storage ls" # e.g. 'gcloud storage' or 'gsutil'
cloud_storage_cp_command="gcloud storage cp" # e.g. 'gcloud storage' or 'gsutil'
cloud_fastq_base_path="gs://mctp-fastq/*" # can include wildcard if cloud storage command can parse them, e.g. 'gs://my-fastq/*'

########################
####   references   ####
########################

ORIG_SGRNA_LIST_FILE=sgrnas.txt

# files derived from ORIG_SGRNA_LIST_FILE will use this base name
SGRNA_LIST_NAME=sgrnas

# only required for bam mode
REFERENCES_DIR=output/references

##################################################################
####   aligner reference flanking sequences (bam mode only)   ####
##################################################################

# these are added to the left and right of reference sgRNA sequences when building the fasta file used to make the bowtie2 index
# current recommendation: leave these values blank

US_SEQ=""    # 5' sequence adjacent to sgRNA, recommended length >= 18 bases
### US_SEQ="taacttgctatttctagctctaaaac"    # 5' sequence adjacent to sgRNA, recommended length >= 18 bases
DS_SEQ=""       # 3' sequence adjacent to sgRNA, recommended length >= 18 bases
### DS_SEQ="cggtgtttcgtcctttccacaag"       # 3' sequence adjacent to sgRNA, recommended length >= 18 bases

#####################################
####   trimming (bam mode only)  ####
#####################################

TRIM_SEQ=TAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC # vector sequence 5' adjacent to the sgRNA sequence
											   # recommended length >= 18 bases, but if too long,
                                               # chances of encountering sequencing errors increases.
                                               # avoid long homopolymers if possible.

#############################################################
####   fastq reorientation (if reorient_fastq is TRUE)   ####
#############################################################

# required if fastq reorientation is required
#   - sequences are searched for in order of appearance in the arrays, R1 first and then R2
#   - if any one sequence is found indicating the read is in the wrong orientation, R1 and R2 are swapped
#   - if no sequence is found, the read will remain in the original orientation
#   - NB: the program does not check for reverse complements of these sequences
#   - at least one sequence is required for r1_seqs; r2_seqs is entirely optional
#   - recommendations
#     - 2 sequences for each r1_seqs and r2_seqs
#     - 16 < sequence length < 24
#     - avoid long homopolymers (higher changes of sequencing errors)
#     - avoid any other regions of the read with high error rates
#     - expected position should be at least several bases from either end of the read

# expected sequences in R1
r1_seqs=(GCTATTTCTAGCTCTAAAAC                 TCCCACTCCTTTCAAGA)
#        ^^^ left adjacent to sgRNA           ^^^ PCR primer (rev) part (5' 5 bases removed)

# expected sequences in R2
r2_seqs=(GTTTTAGAGCTAGAAATAGC                 GAAAGGACGAAACACCG)
#        ^^^ left adjacent to sgRNA (revcomp) ^^^ PCR primer (fwd) part (5' 5 bases removed)

########################################
####   run_pipeline.sh parameters   ####
########################################

# notes:
#  - parameters can be specificed as arguments to run_pipeline.sh as well as defined here
#  - false is recommended for all parameters
#  - cmd line arguments will override these values 

skip_qc=FALSE               # --skip-qc
skip_matrix=FALSE           # --skip-matrix  (skip count/cpm matrix generation; mageck count results will still be generated)
skip_analysis=FALSE         # --skip-analysis (skip comparisons, that is, skip running mageck test)
analysis_only=FALSE         # --analysis-only (skip processing fastq files, generating counts, etc and only run mageck test)
skip_fastq_processing=FALSE # --skip-fastq-processing (if fastq have already been processed)
skip_references=FALSE       # --skip-references (bam mode only; skip reference generation if already generated)

########################
####   QC related   ####
########################

# optional: 
#   - this sequence is used for distance-to-pattern assessment of reads
#   - consider using a sequence closer to the 5 end for speed, possibly one defined in r1_reqs
#   - if not set, then distance-to-pattern will not be calculated

LEN_PATTERN=TCCCACTCCTTTCAAGA
