# this should be sourced

##################################
####   sequencing libraries   ####
##################################

# recommendation: do not modify this section if METADATA_FILE is formatted as described above
# the 3 arrays created below (LIBRARIES, SAMPLE_NAMES, FASTQ_FILES) are required

if [[ "$DOCKER_PATHS" == "TRUE" ]]
then
	FASTQ_DIR=$DOCKER_FASTQ_DIR
fi

if [[ $METADATA_FILE =~ \.ya?ml$ ]]
then

    if ! command -v yq &> /dev/null || ! command -v jq &> /dev/null
    then
        echo "Error: yq and jq are required but not installed." >&2
        exit 1
    fi

	LIBRARIES=($(yq -e '.samples[].library' $METADATA_FILE))
	SAMPLE_NAMES=($(yq -e '.samples[].sample' $METADATA_FILE))
	FASTQ_FILES=($(yq -e '.samples[] | (.fastq_r1 | join(",")) + ";" + (.fastq_r2 | join(","))' $METADATA_FILE))
else
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
fi
