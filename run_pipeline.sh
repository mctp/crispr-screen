#!/bin/bash

set -e  # exit on error

# default values, if not definied in config or via cmd parameters
MODE=fastq
SEARCH_REVCOMP=FALSE
reformat_sgrna_list=FALSE # mostly not used any more, so default is FALSE
skip_references=FALSE
force_references=FALSE
skip_fastq_processing=FALSE
skip_matrix=FALSE
generate_references_cmd_opts=

if [[ -f "config.yml" ]]; then
    if command -v yq &> /dev/null; then
        echo "Parsing config.yml..."
        eval $(yq e 'to_entries | .[] | .key + "=" + (.value | @sh)' config.yml | sed 's/^/export /')
        # yaml variables are all lowercase
        MODE=$mode
        SEARCH_REVCOMP=$search_revcomp
        NCPU=$ncpu
        DOCKER_PATHS=$docker_paths
        WORKING_DIR=$working_dir
        FASTQ_DIR=$fastq_dir
        OUTPUT_DIR=$output_dir
        METADATA_FILE=$metadata_file
        DOCKER_WORKING_DIR=$docker_working_dir
        DOCKER_FASTQ_DIR=$docker_fastq_dir
        ORIG_SGRNA_LIST_FILE=$orig_sgrna_list_file
        SGRNA_LIST_NAME=$sgrna_list_name
        REFERENCES_DIR=$references_dir
        US_SEQ=$us_seq
        DS_SEQ=$ds_seq
        TRIM_SEQ=$trim_seq
        LEN_PATTERN=$len_pattern
    else
        echo "yq is not installed. Falling back to config.sh"
        source config.sh
    fi
else
    source config.sh # fall back to config.sh
fi

source process_metadata.sh

if [[ "$DOCKER_PATHS" != "" ]]
then
    echo "Warning: DOCKER_PATHS is set to $DOCKER_PATHS in config.sh, overriding any autodetection of docker.  This could cause problems."
else
    # detect if we are running in a docker container
    if grep -qE '/docker|/lxc' /proc/1/cgroup 2>/dev/null; then
        DOCKER_PATHS=TRUE
    else
        DOCKER_PATHS=FALSE
    fi
fi

if [[ "$DOCKER_PATHS" == "TRUE" ]]
then
    WORKING_DIR=$DOCKER_WORKING_DIR
    FASTQ_DIR=$DOCKER_FASTQ_DIR
    OUTPUT_DIR=$DOCKER_OUTPUT_DIR
fi

echo "search revcomp: $SEARCH_REVCOMP"

if [[ "$OUTPUT_DIR" == "" || ${#LIBRARIES[@]} -lt 1 ]]
then
    echo "required variables missing. check config.sh and sample_metadata.txt"
    echo "output dir: $OUTPUT_DIR"
    echo "${#SAMPLES[@]} libraries: ${LIBRARIES[@]}"
    exit 1
fi

while [ "$#" -gt 0 ]
do
    case "$1" in
        --skip-matrix )
            skip_matrix=TRUE
            shift 1
            ;;
        --skip-qc )
           skip_qc=TRUE
           shift 1
           ;;
        --skip-analysis )
           skip_analysis=TRUE
           shift 1
           ;;
        --analysis-only )
            analysis_only=TRUE
            shift 1
            ;;
        --skip-fastq-processing )
           skip_fastq_processing=TRUE
           shift 1
           ;;       
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

if [[ "$analysis_only" == "TRUE" ]]
then
    echo "Running analysis section only."
    source mageck_analysis.sh --mode $MODE
    exit
fi

if [[ "$MODE" == "bam" ]]
then
    # check software requirements
    BOWTIE2_PATH=$(which bowtie2)
    CUTADAPT_PATH=$(which cutadapt)
    SAMTOOLS_PATH=$(which samtools)

    if [[  "$SAMTOOLS_PATH" == "" || "$BOWTIE2_PATH" == "" || "$CUTADAPT_PATH" == "" ]]
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
        echo "exiting..."
        exit 1
    fi
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

    if [[ ! "$skip_fastq_processing" == "TRUE" ]]
    then
        # Check if all of the R1 or R2 files exist
        for file in $R1_space $R2_space
        do
            if [[ ! -e "$file" || ! -s "$file" ]]
            then
                echo "Warning: $file does not exist or is empty. Skipping library $LIBRARY."
                # continue to the next library
                continue 2
            fi
        done
    fi

    IN_R1_FASTQ=$OUTPUT_DIR/${LIBRARY}_combined_R1.fq.gz
    IN_R2_FASTQ=$OUTPUT_DIR/${LIBRARY}_combined_R2.fq.gz
    IN_R1_RE_FASTQ=$OUTPUT_DIR/${LIBRARY}_combined_R1.reoriented.fq.gz
    IN_R2_RE_FASTQ=$OUTPUT_DIR/${LIBRARY}_combined_R2.reoriented.fq.gz
    INTERLEAVED_FASTQ=$OUTPUT_DIR/${LIBRARY}_combined_interleaved.fq.gz

    if [[ ! "$skip_fastq_processing" == "TRUE" ]]
    then
        echo "merging $R1_space"
        echo "to: $IN_R1_FASTQ"
        cat $R1_space > $IN_R1_FASTQ

        if [[ "$paired_end_fastq" == "FALSE" ]]
        then
            echo "WARNING: skipping R2 because paired_end_fastq is set to FALSE."
        else
            echo "merging $R2_space"
            echo "to: $IN_R2_FASTQ"
            cat $R2_space > $IN_R2_FASTQ
        fi

        FASTQ_STATS_FILE="$OUTPUT_DIR/fastq_stats.txt"

        if [[ ! -e "$FASTQ_STATS_FILE" ]]
        then
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
            echo "read count:  $NEW_READ_COUNT"
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
            if [[ "$paired_end_fastq" == "FALSE" ]]
            then
                echo "WARNING: skipping reorientation because paired_end_fastq is set to FALSE."
            else

                if [[ ! -e "$IN_R1_RE_FASTQ" || ! -e "$IN_R2_RE_FASTQ" ]]
                then
                    echo "reorienting fastq files (this will take a long time)..."
                    
                    time python reorient_fastq_parallel.py \
                        --sequences-r1 $(IFS=,; echo "${r1_seqs[*]}") \
                        --sequences-r2 $(IFS=,; echo "${r2_seqs[*]}") \
                        --cpus $NCPU \
                        --in-fastq-r1 $IN_R1_FASTQ \
                        --in-fastq-r2 $IN_R2_FASTQ \
                        --out-fastq-r1 $IN_R1_RE_FASTQ \
                        --out-fastq-r2 $IN_R2_RE_FASTQ \
                        --plot \
                        --plot-prefix $OUTPUT_DIR/$SAMPLE

                else
                    echo "reoriented fastq files already exist, skipping reorientation."
                fi
            fi
        fi

    fi

    # QC things
    if [[ "$skip_qc" == "TRUE" ]]
    then
        echo    "WARNING: skipping QC section."    
    else
        echo "--------"
        echo "   R1"
        echo "--------"
        python $WORKING_DIR/find_top_sequences.py \
            --in-fastq $IN_R1_FASTQ \
            --out-fasta $OUTPUT_DIR/${LIBRARY}_R1.fa \
            --out-alignment-clustalw $OUTPUT_DIR/${LIBRARY}_R1.clustalw.aln \
            --sample-size 1000 \
        | tee $OUTPUT_DIR/${LIBRARY}_R1_top_seqs.txt

        if [[ -n "$LEN_PATTERN" ]]; then
            python distance_to_pattern_frequencies.py \
            --in-fastq $IN_R1_FASTQ \
            --pattern $LEN_PATTERN \
            --nreads 10000 \
            --start-at 0 \
            --stop-at 0 \
            --center-sequence \
            --show-threshold 0.005 \
            | tee $OUTPUT_DIR/${LIBRARY}_R1_distance_to_pattern_frequencies.txt
        fi

        if [[ "$paired_end_fastq" == "FALSE" ]]
        then
            echo "WARNING: skipping R2 because paired_end_fastq is set to FALSE."
        else

            echo "--------"
            echo "   R2"
            echo "--------"
            python $WORKING_DIR/find_top_sequences.py \
                --in-fastq $IN_R1_FASTQ \
                --out-fasta $OUTPUT_DIR/${LIBRARY}_R2.fa \
                --out-alignment-clustalw $OUTPUT_DIR/${LIBRARY}_R2.clustalw.aln \
                --sample-size 1000 \
            | tee $OUTPUT_DIR/${LIBRARY}_R2_top_seqs.txt

            if [[ -n "$LEN_PATTERN" ]]; then
                python distance_to_pattern_frequencies.py \
                    --in-fastq $IN_R2_FASTQ \
                    --pattern $LEN_PATTERN \
                    --nreads 10000 \
                    --start-at 0 \
                    --stop-at 0 \
                    --center-sequence \
                    --show-threshold 0.005 \
                | tee $OUTPUT_DIR/${LIBRARY}_R2_distance_to_pattern_frequencies.txt
            fi
        fi

        if [[ -e "$IN_R1_RE_FASTQ" ]]
        then
            echo "-------------------"
            echo "   R1 reoriented"
            echo "-------------------"
            python $WORKING_DIR/find_top_sequences.py \
            --in-fastq $IN_R1_RE_FASTQ \
            --out-fasta $OUTPUT_DIR/${LIBRARY}_R1_reoriented.fa \
            --out-alignment-clustalw $OUTPUT_DIR/${LIBRARY}_R1_reoriented.clustalw.aln \
            --sample-size 1000 \
            | tee $OUTPUT_DIR/${LIBRARY}_R1_reoriented_top_seqs.txt

            if [[ -n "$LEN_PATTERN" ]]; then
                python distance_to_pattern_frequencies.py \
                    --in-fastq $IN_R1_RE_FASTQ \
                    --pattern $LEN_PATTERN \
                    --nreads 10000 \
                    --start-at 0 \
                    --stop-at 0 \
                    --center-sequence \
                    --show-threshold 0.005 \
                | tee $OUTPUT_DIR/${LIBRARY}_R1_reoriented_distance_to_pattern_frequencies.txt
            fi
        else
            echo "WARNING: Reoriented R1 fastq files do not exist. Skipping reoriented QC."
        fi

        if [[ "$paired_end_fastq" == "FALSE" ]]
        then
            echo "WARNING: skipping R2 because paired_end_fastq is set to FALSE."
        else
            if [[ -e "$IN_R2_RE_FASTQ" ]]
            then
                echo "-------------------"
                echo "   R2 reoriented"
                echo "-------------------"
                python $WORKING_DIR/find_top_sequences.py \
                    --in-fastq $IN_R2_RE_FASTQ \
                    --out-fasta $OUTPUT_DIR/${LIBRARY}_R2_reoriented.fa \
                    --out-alignment-clustalw $OUTPUT_DIR/${LIBRARY}_R2_reoriented.clustalw.aln \
                    --sample-size 1000 \
                | tee $OUTPUT_DIR/${LIBRARY}_R2_reoriented_top_seqs.txt

                if [[ -n "$LEN_PATTERN" ]]; then
                    python distance_to_pattern_frequencies.py \
                        --in-fastq $IN_R2_RE_FASTQ \
                        --pattern $LEN_PATTERN \
                        --nreads 10000 \
                        --start-at 0 \
                        --stop-at 0 \
                        --center-sequence \
                        --show-threshold 0.005 \
                    | tee $OUTPUT_DIR/${LIBRARY}_R2_reoriented_distance_to_pattern_frequencies.txt
                fi
            else
                echo "WARNING: Reoriented R2 fastq files do not exist. Skipping reoriented QC."
            fi
        fi
    fi

    if [[ "$reorient_fastq" == "TRUE" && "$MODE" != "bam" && "$paired_end_fastq" != "FALSE" ]]
    then
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
        # check software requirements
        BOWTIE2_PATH=$(which bowtie2)
        CUTADAPT_PATH=$(which cutadapt)
        SAMTOOLS_PATH=$(which samtools)

        if [[  "$SAMTOOLS_PATH" == "" || "$BOWTIE2_PATH" == "" || "$CUTADAPT_PATH" == "" ]]
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
            echo "exiting..."
            exit 1
        fi
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

if [[ "$skip_matrix" == "TRUE" ]]
then
    echo "Skipping matrix generation."
else
    echo
    echo "building count and cpm matrix files..."
    python cpm_matrix.py --config config.sh
    echo
    echo "plotting count distribution histograms..."
    Rscript count_distribution_histograms.R
    echo   
fi

if [[ "$skip_analysis" == "TRUE" ]]
then
    echo "Skipping analysis section."
else
    echo "running mageck analysis..."
    source mageck_analysis.sh --mode $MODE
fi

