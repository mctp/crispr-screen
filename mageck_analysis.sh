#!/bin/bash

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

if [[ "$count_matrix_file" == "" ]]
then
    count_matrix_file=$OUTPUT_DIR/sgrna_count_matrix.txt
fi

echo "experiment name: $experiment_name"
echo "count matrix file: $count_matrix_file"
echo "core essential genes file: $core_essential_genes_file"
echo "non essential genes file: $non_essential_genes_file"
echo "control gene ID: $control_gene_id"
echo "control sample: $control_sample"
echo "control sample min count: $control_sample_min_count"
echo "min sgRNAs per gene: $gene_sgRNA_min"

if [[ "$control_gene_id" == "" ]]
then
    control_genes_file=$non_essential_genes_file
else
    control_genes_file=$OUTPUT_DIR/control_genes.txt
    if [[ ! -f "$control_genes_file" ]]
    then
        echo "$control_gene_id" > $control_genes_file
    fi
fi

if [ -z "${comparisons+x}" ] || [ -z "$comparisons" ]
then
    echo "No comparisons to process."
    return
fi

if [ ! -z "$comparison" ]; then
    comparisons=("$comparisons")
fi



for comparison in "${comparisons[@]}"
do
    echo "comparison: $comparison"
    IFS=':' read -r a b <<< "$comparison"
    # check for valid sample names
    header=$(head -n 1 $count_matrix_file)
    if ! echo "$header" | grep -q -w "$a" || ! echo "$header" | grep -q -w "$b"
    then
        echo "Error: One or both samples ($a, $b) not found in the counts matrix header. Skipping comparison." >&2
        continue
    fi

    if [[ "$normalization_method" == "" ]]
    then
        normalization_method="median"
    elif [[ "$normalization_method" != "median" && "$normalization_method" != "total" && "$normalization_method" != "control" ]]
    then
        echo "Error: Invalid normalization method '$normalization_method'. Valid options are: median, total." >&2
        continue
    fi

    if [[ "$remove_zero_method" == "" ]]
    then
        remove_zero_method="control"
    elif [[ "$remove_zero_method" != "none" && "$remove_zero_method" != "control" && "$remove_zero_method" != "treatment" && "$remove_zero_method" != "both" && "$remove_zero_method" != "any" ]]
    then
        echo "Error: Invalid remove_zero_method '$remove_zero_method'. Valid options are: none, control, treatment, both, any." >&2
        continue
    fi
        
    comparison_name="${a}_vs_${b}"
    echo "comparison name: $comparison_name"
    comparison_dir=$OUTPUT_DIR/$comparison_name
    # prevent overwrite by appending a number to the directory name
    counter=1
    while [ -d "$comparison_dir" ]; do
        comparison_dir="${OUTPUT_DIR}/${comparison_name}_$(printf "%02d" $counter)"
        comparison_prefix=$comparison_dir/$comparison_name
        ((counter++))
    done
    comparison_prefix=$comparison_dir/$comparison_name
    mkdir -p $comparison_dir
    mageck test \
        -k $count_matrix_file \
        -t $a \
        -c $b \
        --control-gene $control_genes_file \
        -n $comparison_prefix \
        --norm-method $normalization_method \
        --remove-zero $remove_zero_method \
        --normcounts-to-file
    summary_file=$comparison_dir/pipeline_variable_summary.txt
    {
        echo "experiment name: $experiment_name"
        echo "count matrix file: $count_matrix_file"
        echo "core essential genes file: $core_essential_genes_file"
        echo "non essential genes file: $non_essential_genes_file"
        echo "control gene ID: $control_gene_id"
        echo "control sample: $control_sample"
        echo "control sample min count: $control_sample_min_count"
        echo "min sgRNAs per gene: $gene_sgRNA_min"
        echo "comparison: $comparison"
        echo "normalization method: $normalization_method"
        echo "remove zero method: $remove_zero_method"
    } > $summary_file
done
