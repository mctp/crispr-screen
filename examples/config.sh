WORKING_DIR=$(pwd)  # path to working dir, this must exist!

OUTPUT_DIR=$WORKING_DIR
METADATA_FILE=$OUTPUT_DIR/sample_metadata.txt # tab delim table with library, sample, fastqs (comma sep)
FASTQ_DIR=  # if metadata gives only file names, do not leave this blank!  if blank, uses metadata verbatim

# references
SGRNA_LIST_NAME=sgrnas_v2
REFERENCES_DIR=$OUTPUT_DIR/references
ORIG_SGRNA_LIST_FILE=$WORKING_DIR/sgrnas_v2.txt

# sequences to add to left and right of sgRNA sequence in the fasta file used to make the bowtie2 index
US_SEQ="taacttgctatttctagctctaaaac"    # 5' sequence adjacent to sgRNA, recommended use at least 18 bases
DS_SEQ="cggtgtttcgtcctttccacaag"       # 3' sequence adjacent to sgRNA, recommended use at least 18 bases
     
    # tagccttattttaacttgctatttctagctctaaaac
# cggtgtttcgtcctttccacaagatat

# trimming
TRIM_SEQ=tagccttattttaacttgctatttctagctctaaaac # trim sequence is 5' vector sequence adjacent to the sgRNA sequence.  it can be long
   #ACGGACTAGCCTTATTTTAACTTGCTATTTCTAGCTCTAAAAC

SEARCH_REVCOMP=TRUE # when generating aligner reference, use reverse complement of the grna sequence (normally leave this FALSE)

# libraries

# filter header and commented lines -- if sample metadata file has no header, create or remove FNR>1 condition!
declare -a LIBRARIES=($(gawk 'BEGIN{FS="\t"; OFS=" "}(FNR>1 && $0!~/^#/){print $1}' $METADATA_FILE))
declare -a SAMPLE_NAMES=($(gawk 'BEGIN{FS="\t"; OFS=" "}(FNR>1 && $0!~/^#/){print $2}' $METADATA_FILE))
declare -a FASTQ_FILES=()
FASTQ_FILES=($(gawk \
	-v fq_path="$FASTQ_DIR" \
		'BEGIN{FS="\t"; OFS=" "
	}(FNR>1 && $0!~/^#/ && $0!~/^[[:space:]]*$/){
		split($3,fq1s,",");
		split($4,fq2s,",");
		if(fq_path != ""){
			new_fq1s = fq_path "/"
			new_fq2s = fq_path "/"
			for(i=1;i<=length(fq1s);i++){
				if(i==1){
					new_fq1s = new_fq1s fq1s[i];
					new_fq2s = new_fq2s fq2s[i];
				}else{
					new_fq1s = new_fq1s "," fq_path "/" fq1s[i];
					new_fq2s = new_fq2s "," fq_path "/" fq2s[i];
				}
			}
		}else{
			new_fq1s = ""
			new_fq2s = ""
			for(i=1;i<=length(fq1s);i++){
				if(i==1){
					new_fq1s = fq1s[i];
					new_fq2s = fq2s[i];
				}else{
					new_fq1s = new_fq1s "," fq1s[i];
					new_fq2s = new_fq2s "," fq2s[i];
				}
			}
		}
		print new_fq1s ";" new_fq2s;
	}' $METADATA_FILE))

# resources
NCPU=12
