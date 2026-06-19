#!/bin/bash

# Extract fastq files from sample metadata (YAML or TXT format)
# Usage: fastq_files_from_sample_metadata.sh [options] [metadata_file]
# Options:
#   --docker    Use yq via docker (default: use local yq if available)
#   --txt [file] Force TXT mode, optionally specify file (default: sample_metadata.txt)

# Default settings
use_docker=false
force_txt=false
txt_file=""
metadata_file=""

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --docker)
            use_docker=true
            shift
            ;;
        --txt)
            force_txt=true
            # Check if next argument exists and doesn't start with --
            if [[ $# -gt 1 && ! "$2" =~ ^-- ]]; then
                txt_file="$2"
                shift 2
            else
                txt_file="sample_metadata.txt"
                shift
            fi
            ;;
        -h|--help)
            echo "Usage: $0 [options] [metadata_file]"
            echo "Options:"
            echo "  --docker      Use yq via docker"
            echo "  --txt [file]  Force TXT mode, optionally specify file"
            echo "  -h, --help    Show this help"
            exit 0
            ;;
        -*)
            echo "Unknown option: $1"
            exit 1
            ;;
        *)
            metadata_file="$1"
            shift
            ;;
    esac
done

# Determine metadata file to use
if [[ "$force_txt" == true ]]; then
    # Force TXT mode
    if [[ "$metadata_file" != "" ]]; then
        # Use provided file if given
        txt_file="$metadata_file"
    fi
    if [[ "$txt_file" == "" ]]; then
        txt_file="sample_metadata.txt"
    fi
    metadata_file="$txt_file"
elif [[ "$metadata_file" == "" ]]; then
    # Auto-detect metadata file
    if [[ -f "sample_metadata.yaml" ]]; then
        metadata_file=sample_metadata.yaml
    elif [[ -f "sample_metadata.txt" ]]; then
        metadata_file=sample_metadata.txt
    else
        echo "Error: Neither sample_metadata.yaml nor sample_metadata.txt found."
        echo "Use --txt to specify a different TXT file or provide a file path."
        exit 1
    fi
fi

if [[ ! -f "$metadata_file" ]]
then
    echo "Error: $metadata_file file not found."
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

# Function to parse YAML manually (basic implementation)
parse_yaml_manually() {
    echo "WARNING: Parsing YAML manually - this may have unintended consequences for complex YAML structures." >&2
    echo "Consider installing yq or using --txt mode for more reliable parsing." >&2
    
    # Basic YAML parsing for the expected structure
    grep -E "^\s*fastq_r[12]:\s*" "$1" | sed 's/^\s*fastq_r[12]:\s*//' | sed 's/^"//' | sed 's/"$//'
}

# Determine file format and processing method
if [[ "$force_txt" == true ]] || [[ "$metadata_file" == *.txt ]]; then
    # Handle TXT format (tab-delimited)
    echo "Using TXT format parsing for: $metadata_file" >&2
    tail -n +2 "$metadata_file" | cut -f2,3 | tr '\t' '\n' | grep -v '^$'
    
elif [[ "$metadata_file" == *.yaml ]] || [[ "$metadata_file" == *.yml ]]; then
    # Handle YAML format
    echo "Using YAML format parsing for: $metadata_file" >&2
    
    # Determine yq method
    if [[ "$use_docker" == true ]]; then
        # User explicitly wants docker
        if docker_yq_works; then
            yq_method="docker"
        elif command_exists yq; then
            echo "WARNING: Docker yq not available. Local yq found - remove --docker to use it" >&2
            yq_method="manual"
        else
            yq_method="manual"
        fi
    else
        # Default: prefer local yq, then docker, then manual
        if command_exists yq; then
            yq_method="local"
        elif docker_yq_works; then
            echo "INFO: Local yq not available. Using docker yq (use --docker to suppress this message)" >&2
            yq_method="docker"
        else
            yq_method="manual"
        fi
    fi
    
    # Execute based on determined method
    case $yq_method in
        "local")
            echo "Using local yq installation" >&2
            yq '.samples[].fastq_r1[], .samples[].fastq_r2[]' "$metadata_file"
            ;;
        "docker")
            echo "Using docker yq" >&2
            docker run --rm -v "$(pwd):/workdir" -w /workdir mikefarah/yq '.samples[].fastq_r1[], .samples[].fastq_r2[]' "$metadata_file"
            ;;
        "manual")
            echo "No yq available - falling back to manual YAML parsing" >&2
            parse_yaml_manually "$metadata_file"
            ;;
    esac
    
else
    echo "Error: Unsupported file format. Use .yaml, .yml, or .txt extension."
    echo "Or use --txt to force TXT mode parsing."
    exit 1
fi


