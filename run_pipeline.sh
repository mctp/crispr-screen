#!/bin/bash

# check software requirements
BOWTIE2_PATH=$(which bowtie2)
CUTADAPT_PATH=$(which cutadapt)
MAGECK_PATH=$(which mageck)
SAMTOOLS_PATH=$(which samtools)

if [[  "$SAMTOOLS_PATH" == "" || "$BOWTIE2_PATH" == "" || "$CUTADAPT_PATH" == "" || "$MAGECK_PATH" == "" ]]
then
	echo "required software is not installed or on the current PATH:"
	if [[  "$SAMTOOLS_PATH" == "" ]]
	then
		echo "samtools"
	fi
	if [[  "$BOWTIE2_PATH" == "" ]]
	then
		echo "bowtie2"
	fi
	if [[  "$CUTADAPT_PATH" == "" ]]
	then
		echo "cutadapt"
	fi
	if [[  "$MAGECK_PATH" == "" ]]
	then
		echo "mageck"
	fi
	echo "exiting..."
	exit 1
fi

source config.sh

if [[ "$WORKING_DIR" == "" || ${#LIBRARIES[@]} -lt 1 ]]
then
	echo "required variables missing. check config.sh and sample_metadata.txt"
	echo "working dir: $WORKING_DIR"
	echo "${#SAMPLES[@]} libraries: ${LIBRARIES[@]}"
	exit 1
fi

reformat_sgrna_list=FALSE # mostly not used any more, so default is FALSE
generate_references_cmd_opts=
skip_references=FALSE

while [ "$#" -gt 0 ]
do
    case "$1" in
        --skip-references )
           skip_references=TRUE
           shift 1
           ;;
        --reformat-sgrna-list )
           reformat_sgrna_list=TRUE
		   generate_references_cmd_opts="$generate_references_cmd_opts --reformat-sgrna-list"
           shift 1
           ;;
        -*|--* )
            echo "unknown option $1" >&2
            exit 1
            ;;
        *)
            args+="$1"
            shift 1
            ;;
    esac
done

# skip ref/index if explicitly set or if output files exist with a test
if [[ "$skip_references" == "TRUE" ]]
then
	echo "skipping reference generataion."
else
	if [[ -e "$REFERENCES_DIR/sgRNAs.fa" ]]
	then
		echo "file exists, skipping reference generation."
	else
		source generate_reference_files.sh $generate_references_cmd_opts
	fi

	if [[ -e "$REFERENCES_DIR/sgRNAs.1.bt2" ]]
	then
		echo "file exists, skipping bt2 index generation."
	else
		source make_bt2_index.sh
	fi 
fi 

for i in ${!LIBRARIES[@]}
do
	LIBRARY="${LIBRARIES[$i]}"
	SAMPLE="${SAMPLE_NAMES[$i]}"

	# change sample name to library id if sample name empty
	if [[ "$SAMPLE" == "" ]]
	then
		SAMPLE="$LIBRARY"
	fi
	IN_R1_FASTQ=$(echo "${FASTQ_FILES[$i]}" | gawk 'BEGIN{FS=","}{print $1}')
	IN_R2_FASTQ=$(echo "${FASTQ_FILES[$i]}" | gawk 'BEGIN{FS=","}{print $2}')

	source trim_reads.sh

	source align_bt2.sh

	source count_mageck.sh

done
