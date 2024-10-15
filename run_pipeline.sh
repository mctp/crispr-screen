#!/bin/bash

set -e  # exit on error

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
        echo "  samtools"
    fi
    if [[  "$BOWTIE2_PATH" == "" ]]
    then
        echo "  bowtie2"
    fi
    if [[  "$CUTADAPT_PATH" == "" ]]
    then
        echo "  cutadapt"
    fi
    if [[  "$MAGECK_PATH" == "" ]]
    then
        echo "  mageck"
    fi
    echo "exiting..."
    exit 1
fi
SEARCH_REVCOMP=FALSE

source config.sh

echo "search revcomp? $SEARCH_REVCOMP"

if [[ "$OUTPUT_DIR" == "" || ${#LIBRARIES[@]} -lt 1 ]]
then
    echo "required variables missing. check config.sh and sample_metadata.txt"
    echo "output dir: $OUTPUT_DIR"
    echo "${#SAMPLES[@]} libraries: ${LIBRARIES[@]}"
    exit 1
fi

reformat_sgrna_list=FALSE # mostly not used any more, so default is FALSE
generate_references_cmd_opts=
skip_references=FALSE
force_references=FALSE

while [ "$#" -gt 0 ]
do
    case "$1" in
        --skip-references )
           skip_references=TRUE
           shift 1
           ;;
        --force-references )
           force_references=TRUE
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

if [[ "$MODE" == "bam" ]]
then
    # skip ref/index if explicitly set or if output files exist with a test
    if [[ "$skip_references" == "TRUE" ]]
    then
        echo "skipping reference generataion."
    else
        if [[ -e "$REFERENCES_DIR/$SGRNA_LIST_NAME.fa"  && "$force_references" == "FALSE" ]]
        then
            echo "file exists, skipping reference generation."
        else
            source generate_reference_files.sh $generate_references_cmd_opts
        fi

        if [[ -e "$REFERENCES_DIR/$SGRNA_LIST_NAME.1.bt2" && "$force_references" == "FALSE" ]]
        then
            echo "file exists, skipping bt2 index generation."
        else
            source make_bt2_index.sh
        fi
    fi
fi


for i in ${!LIBRARIES[@]}
do
    LIBRARY="${LIBRARIES[$i]}"
    SAMPLE="${SAMPLE_NAMES[$i]}"
    echo "$LIBRARY: $SAMPLE"

    # change sample name to library id if sample name empty
    if [[ "$SAMPLE" == "" || FORCE_LIBRARY_AS_NAME == "TRUE" ]]
    then
        SAMPLE="$LIBRARY"
    fi
    IN_R1_FASTQ=$(echo "${FASTQ_FILES[$i]}" | gawk 'BEGIN{FS=";"}{print $1}')
    IN_R2_FASTQ=$(echo "${FASTQ_FILES[$i]}" | gawk 'BEGIN{FS=";"}{print $2}')

    R1_space=$(echo "$IN_R1_FASTQ" | tr ',' ' ')
    R2_space=$(echo "$IN_R2_FASTQ" | tr ',' ' ')

    IN_R1_FASTQ=$OUTPUT_DIR/${LIBRARY}_combined_R1.fq.gz
    IN_R2_FASTQ=$OUTPUT_DIR/${LIBRARY}_combined_R2.fq.gz
    IN_R1_RE_FASTQ=$OUTPUT_DIR/${LIBRARY}_combined_R1.reoriented.fq.gz
    IN_R2_RE_FASTQ=$OUTPUT_DIR/${LIBRARY}_combined_R2.reoriented.fq.gz
    INTERLEAVED_FASTQ=$OUTPUT_DIR/${LIBRARY}_combined_interleaved.fq.gz

    # merge reads, even if only one file
    if [[ ! -e "$IN_R1_FASTQ" ]]
    then
        echo "merging $R1_space"
        echo "to: $IN_R1_FASTQ"
        cat $R1_space > $IN_R1_FASTQ
    fi

    if [[ ! -e "$IN_R2_FASTQ" ]]
    then
        echo "merging $R2_space"
        echo "to: $IN_R2_FASTQ"
        cat $R2_space > $IN_R2_FASTQ
    fi

    FASTQ_STATS_FILE="$OUTPUT_DIR/fastq_stats.txt"

    if [[ ! -e "$FASTQ_STATS_FILE" ]]; then
        echo -e "library\tread_count" > "$FASTQ_STATS_FILE"
    fi

    # Check if the library already has a read count
    READ_COUNT=$(gawk \
        -v lib="$LIBRARY" \
        'BEGIN{
            FS=OFS="\t"
        }{
            if(FNR==1){
                for(i=1; i<=NF; i++){
                    header[$i] = i
                }
            }else{
                if($header["library"] == lib){
                    print $header["read_count"]
                    exit
                }
            }
        }END{
            print -1
        }' "$FASTQ_STATS_FILE")

    if [[ -z "$READ_COUNT" || "$READ_COUNT" -le 0 ]]
    then
        echo "Computing read count for $LIBRARY"
        NEW_READ_COUNT=$(zcat "$IN_R1_FASTQ" | wc -l | gawk '{print $1/4}')
        echo "was:  $READ_COUNT"
        echo "now:  $NEW_READ_COUNT"
        # Report the read count in the fastq_stats file
        if [[ "$READ_COUNT" != "-1" ]]
        then
            # Update existing entry
            gawk \
                -v lib="$LIBRARY" \
                -v count="$NEW_READ_COUNT" \
                'BEGIN{
                    FS=OFS="\t"
                }{
                    if(FNR==1){
                        for(i=1; i<=NF; i++){
                            header[$i] = i
                        }
                    }else{
                        if($header["library"] == lib){
                            $header["read_count"] = count
                        }
                    }
                    print $0
                }' "$FASTQ_STATS_FILE" \
            > "$FASTQ_STATS_FILE.tmp" && mv "$FASTQ_STATS_FILE.tmp" "$FASTQ_STATS_FILE"
        else
            # Add new entry
            gawk \
                -v lib="$LIBRARY" \
                -v count="$NEW_READ_COUNT" \
                'BEGIN{
                    FS=OFS="\t"
                }{
                    if(FNR==1){
                        # find "library" and "read_count" columns
                        for(i=1; i<=NF; i++){
                            header[$i] = i
                        }
                    }
                    print $0
                }END{
                    for(i=1;i<=NF;i++){
                        if(i == header["library"]){
                            if(i==1){
                                printf lib
                            }else{
                                printf OFS lib
                            }
                        }else if(i == header["read_count"]){
                            if(i==1){
                                printf count
                            }else{
                                printf OFS count
                            }
                        }else{
                            if(i==1){
                                printf $i
                            }else{
                                printf OFS $i
                            }
                        }
                    }
                    printf ORS
                }' "$FASTQ_STATS_FILE" \
            > "$FASTQ_STATS_FILE.tmp" && mv "$FASTQ_STATS_FILE.tmp" "$FASTQ_STATS_FILE"
        fi
    fi

    if [[ "$reorient_fastq" == "TRUE" && "$MODE" != "bam" ]]
    then
        if [[ ! -e "$IN_R1_RE_FASTQ" || ! -e "$IN_R2_RE_FASTQ" ]]
        then
            echo "reorienting fastq files (this will take a long time)..."
            time python reorient_fastq_parallel.py \
                --cpus $NCPU \
                --batch-size 1000 \
                "${upstream_seqs[0]}" \
                $IN_R1_FASTQ $IN_R2_FASTQ \
                $IN_R1_RE_FASTQ $IN_R2_RE_FASTQ \
                $INTERLEAVED_FASTQ
        else
            echo "reoriented fastq files already exist, skipping reorientation."
        fi
        # use reoriented fastq files
        IN_R1_FASTQ=$IN_R1_RE_FASTQ
        IN_R2_FASTQ=$IN_R2_RE_FASTQ

    fi

    echo
    echo "r1: $IN_R1_FASTQ"
    echo "r2: $IN_R2_FASTQ"
    echo

    if [[ "$MODE" == "grep" ]]
    then
        # use grep to count sgRNAs in fastq
        source count_grep.sh

    elif [[ "$MODE" == "bam" ]]
    then
        # trim and align reads, then run mageck count with bam input
        source trim_reads.sh
        source align_bt2.sh
        source count_mageck.sh --mode $MODE
    elif [[ "$MODE" == "fastq" ]]
    then
        # run mageck count with fastq input
        source count_mageck.sh --mode $MODE
    fi
done
