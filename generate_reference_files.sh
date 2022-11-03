#!/bin/bash

SGRNA_LIST_PREFIX=$REFERENCES_DIR/$SGRNA_LIST_NAME
SGRNA_LIST_FILE=$SGRNA_LIST_PREFIX.txt
SGRNA_FASTA=$SGRNA_LIST_PREFIX.fa

mkdir -p $REFERENCES_DIR

# convert text to unix style
dos2unix $SGRNA_LIST_FILE

# format for mageck: extract gene name from sgRNA IDs and place into a 3rd column
gawk 'BEGIN{
		FS=OFS="\t"
	}{
		if(FNR==1){
			gene="gene"
		}else{
			gene=gensub(/^sg(.*)_[0-9]+$/,"\\1","g",$1)
		}
		print $1,$2,gene
	}' \
	$ORIG_SGRNA_LIST_FILE \
| sort -k1,1 | uniq \
> $SGRNA_LIST_FILE

# make fasta
gawk 'BEGIN{FS=OFS="\t"}{print ">"$1"\n"$2}' $SGRNA_LIST_FILE > $SGRNA_LIST_NAME.bare_sgrna.fa

# flank grnas with backbone sequence
gawk \
	-v us=$US_SEQ \
	-v ds=$DS_SEQ \
	'{
		if($0!~/^>/){
			$0 = toupper("gttatcaacttgaaaaagtggcaccg") $0 toupper("ctagatcttgagacaaatggc")
		}
		print $0
	}' \
	$SGRNA_LIST_NAME.bare_sgrna.fa \
> $SGRNA_FASTA


