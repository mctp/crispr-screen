#!/bin/bash

SGRNA_LIST_PREFIX=$REFERENCES_DIR/$SGRNA_LIST_NAME
SGRNA_LIST_FILE=$SGRNA_LIST_PREFIX.txt
SGRNA_FASTA=$SGRNA_LIST_PREFIX.fa

mkdir -p $REFERENCES_DIR

echo "search revcomp? $SEARCH_REVCOMP"

# convert text to unix style
dos2unix -n $SGRNA_LIST_FILE $SGRNA_LIST_FILE.tmp
mv $SGRNA_LIST_FILE.tmp $SGRNA_LIST_FILE

reformat_sgrna_list=FALSE # mostly not used any more, so default is FALSE

while [ "$#" -gt 0 ]
do
    case "$1" in
        --reformat-sgrna-list )
           reformat_sgrna_list=TRUE
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

if [[ "$reformat_sgrna_list" == "FALSE" ]]
then
	 
	echo "$ORIG_SGRNA_LIST_FILE"
	source validate_sgrna_list.sh --file $ORIG_SGRNA_LIST_FILE

	cp $ORIG_SGRNA_LIST_FILE $SGRNA_LIST_FILE
else

	source validate_sgrna_list.sh --file $ORIG_SGRNA_LIST_FILE --two-col-format

	#format for mageck: extract gene name from sgRNA IDs and place into a 3rd column
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
fi

# make fa
gawk 'function revcomp(seq) {
		comp["a"]="t"
		comp["c"]="g"
		comp["t"]="a"
		comp["g"]="c"
		split(seq,x,"")
		outstring=""
		for(p=length(x);p>0;p--){
			if(comp[tolower(x[p])]==""){
				outchar=x[p]
			}else{
				outchar=comp[tolower(x[p])]
				if(tolower(x[p]) != x[p]){ # was uppercase
					outchar=toupper(outchar)
				}
			}
			outstring = outstring outchar;
		}
		return(outstring)
	}
	BEGIN{
		FS=OFS="\t"
	}(FNR>1){
		rc = revcomp($2)
		printf ">" $1 ORS rc ORS
	}' $SGRNA_LIST_FILE \
> $REFERENCES_DIR/$SGRNA_LIST_NAME.bare_sgrna.rc.fa

gawk 'BEGIN{FS=OFS="\t"}(FNR>1){print ">"$1"\n"$2}' $SGRNA_LIST_FILE > $REFERENCES_DIR/$SGRNA_LIST_NAME.bare_sgrna.fa

if [[ "$SEARCH_REVCOMP" == "TRUE" ]]
then
	BARE_SGRNA_FA=$REFERENCES_DIR/$SGRNA_LIST_NAME.bare_sgrna.rc.fa
else
	BARE_SGRNA_FA=$REFERENCES_DIR/$SGRNA_LIST_NAME.bare_sgrna.fa
fi
# flank grnas with backbone sequence
gawk \
	-v us=$US_SEQ \
	-v ds=$DS_SEQ \
	'{
		if($0!~/^>/){
			$0 = toupper(us) $0 toupper(ds)
		}
		print $0
	}' \
	$BARE_SGRNA_FA \
> $SGRNA_FASTA


