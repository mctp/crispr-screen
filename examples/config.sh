## configuration file for Mahnoor Gondal tutorial ##

## resources section ##
NCPU=32

## paths section ##
# NOTE: when running in docker, these are relative to docker container, not host

WORKING_DIR=/repo # path to code base (git repo)

FASTQ_DIR=/input # put fastq files in FASTQ_DIR
                 # if metadata provides full path to fastq dir, this can be left blank

METADATA_FILE=$WORKING_DIR/sample_metadata.txt # tab delim table with 4 columns: library, sample, fastq_r1, fastrq_r2
                                               # header is required; comments (lines beginning with #) are ignored
                                               # NOTE: where R1 and R2 have multiple files use comma-separated list

OUTPUT_DIR=/output # this dir must exist. if using docker, write permissions may need to be explicitly set

## references section ##

REFERENCES_DIR=$OUTPUT_DIR/references
ORIG_SGRNA_LIST_FILE=$WORKING_DIR/sgRNAs.txt
SGRNA_LIST_NAME=sgRNAs

# aligner reference flanking sequences
# these are added to the left and right of reference sgRNA sequences when building the fasta file used to make the bowtie2 index
US_SEQ="taacttgctatttctagctctaaaac"    # 5' sequence adjacent to sgRNA, recommended length >= 18 bases
DS_SEQ="cggtgtttcgtcctttccacaag"       # 3' sequence adjacent to sgRNA, recommended length >= 18 bases

## trimming section ##

TRIM_SEQ=tagccttattttaacttgctatttctagctctaaaac # trim sequence is 5' vector sequence adjacent to the sgRNA sequence
                                               # recommended length > 18 bases, but it can be longer

SEARCH_REVCOMP=TRUE # when generating aligner reference, use reverse complement of the sgrna sequence
                    # try changing this if alignment rate is low

## libraries section ##

# recommendation is to not modify this section if METADATA_FILE is formatted as described above
# the 3 arrays below are required

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

