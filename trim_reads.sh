#!/bin/bash

SAMPLE_DIR=$OUTPUT_DIR/${SAMPLE}_bam
SAMPLE_PREFIX=$SAMPLE_DIR/${SAMPLE}_bam
OUT_R1_FASTQ=$SAMPLE_PREFIX.R1.fq
OUT_R2_FASTQ=$SAMPLE_PREFIX.R2.fq
OUT_R1_UNTRIMMED_FASTQ=$SAMPLE_PREFIX.R1.untrimmed.fq
OUT_R2_UNTRIMMED_FASTQ=$SAMPLE_PREFIX.R2.untrimmed.fq
OUT_COMBINED_FASTQ=$SAMPLE_PREFIX.fq

mkdir -p $SAMPLE_DIR

if [[ "$(which parsort 2>/dev/null)" != "" ]]
then
	SORT="parsort"
else
	SORT="sort"
fi

echo "using $SORT to sort."

# run cutadapt for each R1/R2 then combine
# cutadapt -m 18 -j $NCPU -O 18 -g $TRIM_SEQ -l 20 --discard-untrimmed -o $OUT_R1_FASTQ $IN_R1_FASTQ
# cutadapt -m 18 -j $NCPU -O 18 -g $TRIM_SEQ -l 20 --discard-untrimmed -o $OUT_R2_FASTQ $IN_R2_FASTQ

cutadapt -m 18 -j $NCPU -O 18 -g $TRIM_SEQ -l 20 --untrimmed-output $OUT_R1_UNTRIMMED_FASTQ -o $OUT_R1_FASTQ $IN_R1_FASTQ
cutadapt -m 18 -j $NCPU -O 18 -g $TRIM_SEQ -l 20 --untrimmed-output $OUT_R2_UNTRIMMED_FASTQ -o $OUT_R2_FASTQ $IN_R2_FASTQ
cat $OUT_R1_FASTQ $OUT_R2_FASTQ > $OUT_COMBINED_FASTQ


# check that no read names are duplicated --- no read pair should BOTH pass trimming filter!

echo "checking for duplicated read names in output fastq...."
N_DUP_NAME=$(gawk '(FNR % 4 == 1 ){print $1}' $OUT_COMBINED_FASTQ | $SORT | uniq -c | sed -E 's/^[ \t]+//g;s/[ \t]+/\t/g' | gawk 'BEGIN{FS=OFS="\t"; ct=0}{if($1>1){ct++}}END{print ct}')
echo  "$N_DUP_NAME duplicated read names.  WARNING: if this number is > 0, some DNA fragments may be double-counted!"

#notes:
#   -m 18 set along with -l 20 and --discard-untrimmed results in all reads being 18-20 base only
#   -O 18 means 
