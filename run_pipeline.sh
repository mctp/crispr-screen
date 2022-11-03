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

# source generate_reference_files.sh

# source make_bt2_index.sh

#for i in ${!LIBRARIES[@]}
for i in 0
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

#	source trim_reads.sh

#	source align_bt2.sh

	source count_mageck.sh

done
