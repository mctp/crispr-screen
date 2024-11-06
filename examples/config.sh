#------------------------------------------------#
##################################################
####   CRISPR screen pipeline configuration   ####
##################################################
#------------------------------------------------#

MODE=fastq # fastq (preferred), bam, grep (experimentation only)

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

METADATA_FILE=$WORKING_DIR/sample_metadata.txt # see sample_metadata.txt and sample_metadata.yaml in examples/
                                               # for valid metadata files


#########################
####   cloud paths   ####
#########################

cloud_storage_command="gcloud storage" # e.g. 'gcloud storage' or 'gsutil'
cloud_storage_ls_command="$cloud_storage_command ls" # e.g. 'gcloud storage' or 'gsutil'
cloud_storage_cp_command="$cloud_storage_command cp" # e.g. 'gcloud storage' or 'gsutil'
cloud_fastq_base_path="gs://mctp-fastq/*" # can include wildcard if cloud storage command can parse them, e.g. 'gs://my-fastq/*'

########################
####   references   ####
########################

ORIG_SGRNA_LIST_FILE=$WORKING_DIR/sgrnas.txt

# files derived from ORIG_SGRNA_LIST_FILE will use this base name
SGRNA_LIST_NAME=sgrnas

# only required for bam mode
REFERENCES_DIR=$OUTPUT_DIR/references

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

TRIM_SEQ=TAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC # TRIM SEQUENCE IS 5' VECTOR SEQUENCE ADJACENT TO THE SGRNA SEQUENCE
											   # RECOMMENDED LENGTH > 18 BASES, BUT IT CAN BE LONGER

#############################################################
####   fastq reorientation (if reorient_fastq is TRUE)   ####
#############################################################

# sequences are searched for in order of appearance in the arrays, R1 first and then R2
# recommendation is to select sequence adjacent to sgRNA and/or PCR primers
# when using PCR primers, leave off first 3-4 bases on the ends of the amplicon
# the code will not check reverse complement of these sequences!

# expected sequences in R1
r1_seqs=(GCTATTTCTAGCTCTAAAAC                 TCCCACTCCTTTCAAGA)
#        ^^^ left adjacent to sgRNA           ^^^ PCR primer (rev) part (5' 5 bases removed)

# expected sequences in R2
r2_seqs=(GTTTTAGAGCTAGAAATAGC                 GAAAGGACGAAACACCG)
#        ^^^ left adjacent to sgRNA (revcomp) ^^^ PCR primer (fwd) part (5' 5 bases removed)

########################
####   QC related   ####
########################

# this sequence is used for distance-to-pattern assessment of reads
LEN_PATTERN=AGTTACGCCAAGC
