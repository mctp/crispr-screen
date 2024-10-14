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
#trim5_lengths=120 # comma-separated numbers or leave empty for default (recommended)
#keep_tmp=TRUE  # FALSE recommended
#save_unmapped=TRUE # FALSE recommended

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

# NOTE: when running in docker, these are relative to docker container, not host

WORKING_DIR=/repo # path to code base (git repo)

FASTQ_DIR=/input # put fastq files in FASTQ_DIR
                 # if metadata provides full path to fastq dir, this can be left blank

METADATA_FILE=$WORKING_DIR/sample_metadata.txt # tab delim table with 4 columns: library, sample, fastq_r1, fastrq_r2
                                               # header is required; comments (lines beginning with #) are ignored
                                               # NOTE: where R1 and R2 have multiple files use comma-separated list

OUTPUT_DIR=/output # this dir must exist. if using docker, write permissions may need to be explicitly set

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
# recommendation is to leave these values blank

US_SEQ=""    # 5' sequence adjacent to sgRNA, recommended length >= 18 bases
### US_SEQ="taacttgctatttctagctctaaaac"    # 5' sequence adjacent to sgRNA, recommended length >= 18 bases
DS_SEQ=""       # 3' sequence adjacent to sgRNA, recommended length >= 18 bases
### DS_SEQ="cggtgtttcgtcctttccacaag"       # 3' sequence adjacent to sgRNA, recommended length >= 18 bases

#####################################
####   trimming (bam mode only)  ####
#####################################

TRIM_SEQ=tagccttattttaacttgctatttctagctctaaaac # trim sequence is 5' vector sequence adjacent to the sgRNA sequence
                                               # recommended length > 18 bases, but it can be longer

#############################################################
####   fastq reorientation (if reorient_fastq is TRUE)   ####
#############################################################

# sequences are searched for in order of appearance in the arrays, R1 first and then R2
# recommendation is to select sequence adjacent to sgRNA and/or PCR primers 
# when using PCR primers, leave off first 3-4 bases on the ends of the amplicon
# the code will not check reverse complement of these sequences!

# upstream sequences define new R1
upstream_seqs=(cttgctatttctagctctaaaac)
# downstream sequences define new R2
downstream_seqs=(cggtgtttcgtcctttccacaag)

##################################
####   sequencing libraries   ####
##################################

# recommendation: do not modify this section if METADATA_FILE is formatted as described above
# the 3 arrays created below are required

LIBRARIES=($(gawk 'BEGIN{FS="\t"; OFS=" "}(FNR>1 && $0!~/^#/){print $1}' $METADATA_FILE))
SAMPLE_NAMES=($(gawk 'BEGIN{FS="\t"; OFS=" "}(FNR>1 && $0!~/^#/){print $2}' $METADATA_FILE))
# fastq format is <R1_file1>,<R1_file2>;<R2_file1>,<R2_file2>
FASTQ_FILES=($(gawk \
	-v fq_path="$FASTQ_DIR" \
		'BEGIN{FS="\t"; OFS=" "
	}(FNR>1 && $0!~/^#/ && $0!~/^[[:space:]]*$/){
		split($3,fq1s,",");
		split($4,fq2s,",");
		if(fq_path != ""){
			fq_path_slash = fq_path "/"
		}else{
			fq_path_slash = ""
		}
		for(i=1;i<=length(fq1s);i++){
			if(i==1){
				new_fq1s = fq_path_slash fq1s[i];
				new_fq2s = fq_path_slash fq2s[i];
			}else{
				new_fq1s = new_fq1s "," fq_path_slash fq1s[i];
				new_fq2s = new_fq2s "," fq_path_slash fq2s[i];
			}
		}
		print new_fq1s ";" new_fq2s;
	}' $METADATA_FILE))

