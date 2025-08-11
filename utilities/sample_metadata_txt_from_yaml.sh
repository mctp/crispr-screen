#!/bin/bash

# Convert YAML sample metadata to TXT format
# Usage: sample_metadata_txt_from_yaml.sh [options] [yaml_file]
# Options:
#   --docker    Use yq via docker (default: use local yq if available)

# Default settings
use_docker=false
yaml_file=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --docker)
            use_docker=true
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [options] [yaml_file]"
            echo "Options:"
            echo "  --docker      Use yq via docker"
            echo "  -h, --help    Show this help"
            exit 0
            ;;
        -*)
            echo "Unknown option: $1"
            exit 1
            ;;
        *)
            yaml_file="$1"
            shift
            ;;
    esac
done

if [[ "$yaml_file" == "" ]]
then
    yaml_file=sample_metadata.yaml
fi

if [[ ! -f "$yaml_file" ]]
then
    echo "Error: $yaml_file file not found."
    exit 1
fi

# Function to check if command exists
command_exists() {
    command -v "$1" >/dev/null 2>&1
}

# Function to check if docker yq works
docker_yq_works() {
    docker run --rm mikefarah/yq --version >/dev/null 2>&1
}

# Determine yq method
if [[ "$use_docker" == true ]]; then
    # User explicitly wants docker
    if docker_yq_works; then
        yq_method="docker"
        echo "Using docker yq" >&2
    elif command_exists yq; then
        echo "WARNING: Docker yq not available. Local yq found - remove --docker to use it" >&2
        echo "Falling back to local yq" >&2
        yq_method="local"
    else
        echo "ERROR: Neither docker yq nor local yq available" >&2
        exit 1
    fi
else
    # Default: prefer local yq, then docker
    if command_exists yq; then
        yq_method="local"
        echo "Using local yq installation" >&2
    elif docker_yq_works; then
        echo "INFO: Local yq not available. Using docker yq (use --docker to suppress this message)" >&2
        yq_method="docker"
    else
        echo "ERROR: Neither local yq nor docker yq available" >&2
        echo "Please install yq or ensure docker is available with mikefarah/yq image" >&2
        exit 1
    fi
fi

# Write header
printf "library\tsample\tfastq_r1\tfastq_r2\n" > sample_metadata.txt

# Execute yq command based on determined method
tab=$'\t'
case $yq_method in
    "local")
        if ! yq -r ".samples[] | [.library, .sample, (.fastq_r1 | join(\",\")), (.fastq_r2 | join(\",\"))] | join(\"$tab\")" "$yaml_file" >> sample_metadata.txt; then
            echo "Error: Failed to process YAML with local yq"
            exit 1
        fi
        ;;
    "docker")
        if ! docker run --rm -v "$(pwd):/workdir" -w /workdir mikefarah/yq -r ".samples[] | [.library, .sample, (.fastq_r1 | join(\",\")), (.fastq_r2 | join(\",\"))] | join(\"$tab\")" "$yaml_file" >> sample_metadata.txt; then
            echo "Error: Failed to process YAML with docker yq"
            exit 1
        fi
        ;;
esac

# Check if the file was created successfully and has content
if [[ ! -f "sample_metadata.txt" ]]; then
    echo "Error: sample_metadata.txt was not created."
    exit 1
elif [[ ! -s "sample_metadata.txt" ]]; then
    echo "Error: sample_metadata.txt is empty."
    exit 1
else
    echo "sample_metadata.txt created successfully."
fi  

# Validate the output file
# Check that first line (header) has 4 columns
header_cols=$(head -n 1 sample_metadata.txt | awk -F'\t' '{print NF}')
if [[ $header_cols -ne 4 ]]; then
    echo "Error: Header line has $header_cols columns, expected 4."
    exit 1
fi

# Check that all data lines have 4 columns  
data_line_count=$(tail -n +2 sample_metadata.txt | wc -l)
if [[ $data_line_count -eq 0 ]]; then
    echo "Error: No data lines found in sample_metadata.txt (only header present)."
    exit 1
fi

# Check column count consistency in data lines
inconsistent_lines=$(tail -n +2 sample_metadata.txt | awk -F'\t' '{print NF}' | sort -u | wc -l)
if [[ $inconsistent_lines -ne 1 ]]; then
    echo "Error: Data lines have inconsistent number of columns."
    tail -n +2 sample_metadata.txt | awk -F'\t' '{print NR+1 ": " NF " columns"}' | head -5
    exit 1
fi

# Check that all data lines have exactly 4 columns
data_cols=$(tail -n +2 sample_metadata.txt | awk -F'\t' '{print NF}' | sort -u)
if [[ $data_cols -ne 4 ]]; then
    echo "Error: Data lines have $data_cols columns, expected 4."
    exit 1
fi

echo "sample_metadata.txt is valid ($data_line_count data rows)."
