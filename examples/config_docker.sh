# for docker use, dirs are relative to the docker container, not the host!

WORKING_DIR=/repo  # path to working dir, this must exist!

METADATA_FILE=$WORKING_DIR/sample_metadata.txt # tab delim table with library, sample, fastqs (comma sep)

OUTPUT_DIR=/output  
FASTQ_DIR=/input   # if metadata gives only file names, do not leave this blank!  if blank, uses metadata verbatim

# references
SGRNA_LIST_NAME=sgRNAs
REFERENCES_DIR=$OUTPUT_DIR/references
ORIG_SGRNA_LIST_FILE=$WORKING_DIR/human.sgRNAs.txt

# sequences to add to left and right of sgRNA sequence in the fasta file used to make the bowtie2 index
US_SEQ="gttatcaacttgaaaaagtggcaccg"  # 3' sequence adjacent to sgRNA, recommended use at least 18 bases
DS_SEQ="ctagatcttgagacaaatggc"       # 5' sequence adjacent to sgRNA, recommended use at least 18 bases

# trimming
TRIM_SEQ=tatatcttgtggaaaggacgaaacaccg # trim sequence is 5' vector sequence adjacent to the sgRNA sequence.  it can be long

# libraries

# filter header and commented lines -- if sample metadata file has no header, create or remove FNR>1 condition!
LIBRARIES=($(gawk 'BEGIN{FS="\t"; OFS=" "}(FNR>1 && $0!~/^#/){print $1}' $METADATA_FILE))
SAMPLE_NAMES=($(gawk 'BEGIN{FS="\t"; OFS=" "}(FNR>1 && $0!~/^#/){print $2}' $METADATA_FILE))
FASTQ_FILES=($(gawk \
	-v fq_path="$FASTQ_DIR" \
		'BEGIN{FS="\t"; OFS=" "
	}(FNR>1 && $0!~/^#/){
		if(fq_path!=""){
			split($3,fqs,",");
			new_fq = fq_path "/" fqs[1] "," fq_path "/" fqs[2];
			print new_fq
		}else{
			print $3;
		}
	}' $METADATA_FILE))

# resources
NCPU=12
