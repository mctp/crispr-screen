#!/bin/bash

echo "ignore R2? $ignore_r2"

while [ "$#" -gt 0 ]
do
    case "$1" in
        --mode )
           mageck_mode=$2
           shift 2
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

if [[ "$mageck_mode" == "bam" || "$mageck_mode" == "" ]]
then
    SAMPLE_DIR=$OUTPUT_DIR/${SAMPLE}
else
    SAMPLE_DIR=$OUTPUT_DIR/${SAMPLE}_${mageck_mode}
fi

SAMPLE_PREFIX=$SAMPLE_DIR/$SAMPLE
OUT_BAM=$SAMPLE_PREFIX.aligned.bam
MAGECK_OUT_PREFIX=$SAMPLE_PREFIX

mkdir -p $SAMPLE_DIR


if [[ "$mageck_mode" == "bam" || "$mageck_mode" == "" ]]
then
    SGRNA_LIST_PREFIX=$REFERENCES_DIR/$SGRNA_LIST_NAME # also same as bt2 index prefix
    SGRNA_LIST_FILE=$SGRNA_LIST_PREFIX.txt
    SGRNA_FASTA=$SGRNA_LIST_PREFIX.fa
    mkdir -p $REFERENCES_DIR
    # count reads for sgRNAs after alignment
    mageck count -l $SGRNA_LIST_FILE -n $MAGECK_OUT_PREFIX --sample-label "$SAMPLE" --fastq $OUT_BAM
elif [[ "$mageck_mode" == "fastq" ]]
then
    SGRNA_LIST_FILE=$SGRNA_LIST_NAME.txt
    if [[ "$ignore_r2" == "TRUE" ]]
    then
        fastq2_opt=""
        count_pair_opt=""
        other_opt="--trim-5 106 --sgrna-len 20 --unmapped-to-file --keep-tmp" 
    else
        fastq2_opt="--fastq-2 $IN_R2_FASTQ"
        count_pair_opt="--count-pair True"
        other_opt=""
    fi

    if [[ "$SEARCH_REVCOMP" == "TRUE" ]]
    then
        revcomp_opt="--reverse-complement"
    fi
    # count reads for sgRNAs from fastq
    mageck count \
        -l "$SGRNA_LIST_FILE" \
        -n "$MAGECK_OUT_PREFIX" \
        --sample-label "$SAMPLE" \
        $revcomp_opt \
        --fastq "$IN_R1_FASTQ" \
        $fastq2_opt \
        $count_pair_opt \
        $other_opt
fi
