# source this code


function revcomp {
    # input string of actg, case insensitive, return reversed complement
    gawk \
        -v input="$1" \
        'BEGIN{
            comp["a"]="t"
            comp["c"]="g"
            comp["t"]="a"
            comp["g"]="c"
            split(input,x,"")
            for(p=length(x);p>0;p--){
                if(comp[tolower(x[p])]==""){
                    outchar=x[p]
                }else{
                    outchar=comp[tolower(x[p])]
                    if(tolower(x[p]) != x[p]){ # was uppercase
                        outchar=toupper(outchar)
                    }
                }
                printf outchar
            }
            printf ORS
        }'
}

IN_R1_FASTQ_UNZIPPED="${IN_R1_FASTQ%.gz}"
IN_R2_FASTQ_UNZIPPED="${IN_R2_FASTQ%.gz}"
IN_R1_FASTQ_SEQ="${IN_R1_FASTQ%.fq.gz}.txt"
IN_R2_FASTQ_SEQ="${IN_R2_FASTQ%.fq.gz}.txt"

if [[ -e "$OUTPUT_DIR/${SAMPLE}_${MODE}/${SAMPLE}_count.txt" ]]
then
    rm $OUTPUT_DIR/${SAMPLE}_${MODE}/${SAMPLE}_count.txt
fi

# if [[ ! -e "$IN_R1_FASTQ_UNZIPPED" ]]
# then
#     echo "unzipping $IN_R1_FASTQ..."
#     gunzip -c $IN_R1_FASTQ > $IN_R1_FASTQ_UNZIPPED
# fi

if [[ ! -e "$IN_R1_FASTQ_SEQ" ]]
then
    echo "R1: unzipping fastq and keeping sequences..."
    gunzip -c $IN_R1_FASTQ | gawk 'FNR % 4 == 2' > $IN_R1_FASTQ_SEQ
fi

# if [[ ! -e "$IN_R2_FASTQ_UNZIPPED" ]]
# then
#     echo "unzipping $IN_R2_FASTQ..."
#     gunzip -c $IN_R2_FASTQ > $IN_R2_FASTQ_UNZIPPED
# fi

if [[ ! -e "$IN_R2_FASTQ_SEQ" ]]
then
    echo "R2: unzipping fastq and keeping sequences..."
    gunzip -c $IN_R2_FASTQ | gawk 'FNR % 4 == 2' > $IN_R2_FASTQ_SEQ
fi

echo "sample: $SAMPLE"
SGRNA_LIST_FILE=$SGRNA_LIST_NAME.txt
echo "counting sgRNAs in $SGRNA_LIST_FILE..."
sgrna_ids=($(gawk 'BEGIN{FS=OFS="\t"}(FNR>1){print $1}' $SGRNA_LIST_FILE))
sgrna_seqs=($(gawk 'BEGIN{FS=OFS="\t"}(FNR>1){print $2}' $SGRNA_LIST_FILE))
sgrna_targets=($(gawk 'BEGIN{FS=OFS="\t"}(FNR>1){print $3}' $SGRNA_LIST_FILE))

mkdir -p $OUTPUT_DIR/${SAMPLE}_${MODE}
echo -e "sgrna_id\tsgrna_sequence\tsgrna_target\tR1_count\tR2_count\tR1_rc_count\tR2_rc_count" > $OUTPUT_DIR/${SAMPLE}_${MODE}/${SAMPLE}_count.txt



function count_sgrna() {
    local sgrna_id=$1
    local sgrna_seq=$2
    local sgrna_target=$3

    # echo "  counting $sgrna_id ($sgrna_seq), targeting $sgrna_target ..."
    r1_count=$(grep -c "$sgrna_seq" $IN_R1_FASTQ_SEQ)
    r2_count=$(grep -c "$sgrna_seq" $IN_R2_FASTQ_SEQ)
    r1_rc_count=$(grep -c $(revcomp $sgrna_seq) $IN_R1_FASTQ_SEQ)
    r2_rc_count=$(grep -c $(revcomp $sgrna_seq) $IN_R2_FASTQ_SEQ)
    echo -e "$sgrna_id\t$sgrna_seq\t$sgrna_target\t$r1_count\t$r2_count\t$r1_rc_count\t$r2_rc_count" >> $OUTPUT_DIR/${SAMPLE}_${MODE}/${SAMPLE}_count.txt
}

# Combine the arrays into a single array of tuples
combined_array=()
for i in "${!sgrna_ids[@]}"; do
    combined_array+=("${sgrna_ids[$i]},${sgrna_seqs[$i]},${sgrna_targets[$i]}")
done

# exports required for parallel
export -f count_sgrna
export -f revcomp
export IN_R1_FASTQ_SEQ
export IN_R2_FASTQ_SEQ
export OUTPUT_DIR
export SAMPLE
export MODE

parallel -j $NCPU --colsep ',' count_sgrna ::: "${combined_array[@]}"
# parallel -j $NCPU count_sgrna ::: "${sgrna_ids[@]}" ::: "${sgrna_seqs[@]}" ::: "${sgrna_targets[@]}"
