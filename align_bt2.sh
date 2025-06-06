#!/bin/bash

# output names
SGRNA_LIST_PREFIX=$REFERENCES_DIR/$SGRNA_LIST_NAME # also same as bt2 index prefix
SGRNA_FASTA=$SGRNA_LIST_PREFIX.fa

SAMPLE_DIR=$OUTPUT_DIR/${SAMPLE}_bam
SAMPLE_PREFIX=$SAMPLE_DIR/${SAMPLE}_bam
OUT_BAM=$SAMPLE_PREFIX.aligned.bam
OUT_STATS=$SAMPLE_PREFIX.aligned.stats

OUT_R1_FASTQ=$SAMPLE_PREFIX.R1.fq
OUT_R2_FASTQ=$SAMPLE_PREFIX.R2.fq
OUT_COMBINED_FASTQ=$SAMPLE_PREFIX.fq
OUT_UNALIGNED_FASTQ_GZ=$SAMPLE_PREFIX.unaligned.fq.gz

mkdir -p $SAMPLE_DIR

echo "begin bowtie2 alignment..."

# align
bowtie2 -p $NCPU -x $SGRNA_LIST_PREFIX -U $OUT_COMBINED_FASTQ --norc --un-gz $OUT_UNALIGNED_FASTQ_GZ | samtools view -@ $NCPU -bS - > $OUT_BAM
#bowtie2 -p $NCPU -x $SGRNA_LIST_PREFIX -U $OUT_COMBINED_FASTQ | samtools view -@ $NCPU -bS - > $OUT_BAM

# stats
samtools stats $OUT_BAM > $OUT_STATS
