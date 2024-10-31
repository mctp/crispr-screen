#!/bin/bash

# Parse command line arguments
skip_file_check=FALSE
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --skip-file-check) skip_file_check=TRUE ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# Source the configuration file
source config.sh

echo "downloading using: $cloud_storage_command"
echo "from base path: $cloud_fastq_base_path"

# Process the FASTQ_FILES array
for file_list in "${FASTQ_FILES[@]}"
do
    # Split the list by semicolons and commas
    IFS=';, ' read -r -a files <<< "$file_list"
    
    for file in "${files[@]}"
    do
        # Extract the basename of the file
        base_file=$(basename "$file")
        
        if [[ "$skip_file_check" == "TRUE" ]]
        then
            if [[ "$base_file" =~ mctp_SI_[0-9]+_([^_]+)_[0-9]+_[0-9]+\.fq\.gz ]]
            then
                flowcell="${BASH_REMATCH[1]}"
                result="${cloud_fastq_base_path/\*/$flowcell}/$base_file"
            else
                result="$cloud_fastq_base_path/$base_file"
            fi
            echo "* $result"
        else
            if [[ "$base_file" =~ mctp_SI_[0-9]+_([^_]+)_[0-9]+_[0-9]+\.fq\.gz ]]
            then
            flowcell="${BASH_REMATCH[1]}"
            result=$($cloud_storage_ls_command "${cloud_fastq_base_path/\*/$flowcell}/$base_file")
            else
            result=$($cloud_storage_ls_command "$cloud_fastq_base_path/$base_file")
            fi

            # Check if the result contains more than one file path
            if [[ "$result" == *$'\n'* ]]
            then
            echo "Warning: Multiple files found for $base_file. Listing all files:"
            files_found=($(echo "$result"))
            for i in "${!files_found[@]}"
            do
                if [[ $i -eq 0 ]]
                then
                echo "* ${files_found[$i]}"
                else
                echo "  ${files_found[$i]}"
                fi
            done
            result="${files_found[0]}"
            else
            echo "* $result"
            fi
        fi
        if [[ "$result" != "" ]]
        then
            # Download the file to the specified directory
            $cloud_storage_cp_command "$result" "$FASTQ_DIR"
        else
            echo "File $base_file does not exist in the bucket."
        fi
    done
done