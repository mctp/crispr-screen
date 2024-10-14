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

if [[ "$mageck_mode" == "" ]]
then
    SAMPLE_DIR=$OUTPUT_DIR/${SAMPLE}
    SAMPLE_PREFIX=$SAMPLE_DIR/$SAMPLE
else
    SAMPLE_DIR=$OUTPUT_DIR/${SAMPLE}_${mageck_mode}
    SAMPLE_PREFIX=$SAMPLE_DIR/${SAMPLE}_${mageck_mode}
fi

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

    if [[ "$trim5_lengths" == "" ]]
    then
        other_opts=""
    else
        other_opts="--trim-5 $trim5_lengths"
    fi

    if [[ "$keep_tmp" == "TRUE" ]]
    then
        other_opts="$other_opts --keep-tmp"
    fi

    if [[ "$save_unmapped" == "TRUE" ]]
    then
        other_opts="$other_opts --unmapped-to-file"
    fi

    if [[ "$trim5_lengths" != "" || "$keep_tmp" == "TRUE" || "$save_unmapped" == "TRUE" ]]
    then
        echo "forcing sgrna length to 20 because options trim5_lengths, keep_tmp, and/or save_unmapped are set in config."
        other_opts="$other_opts --sgrna-len 20"
    fi

    if [[ "$ignore_r2" != "TRUE" ]]
    then
        other_opts="$other_opts --fastq-2 $IN_R2_FASTQ"
        other_opts="$other_opts --count-pair True"
    fi

    if [[ "$SEARCH_REVCOMP" == "TRUE" ]]
    then
        other_opts="$other_opts --reverse-complement"
    fi
    # count reads for sgRNAs from fastq
    mageck count \
        -l "$SGRNA_LIST_FILE" \
        -n "$MAGECK_OUT_PREFIX" \
        --sample-label "$SAMPLE" \
        --fastq "$IN_R1_FASTQ" \
        $other_opts
fi
