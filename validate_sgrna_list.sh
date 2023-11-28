# this script should be sourced by other scripts

# purpose: validate the file format of the input sgrna list file

two_col_format=FALSE
file=

while [ "$#" -gt 0 ]
do
    case "$1" in
        --two-col-format )
           two_col_format=TRUE
           shift 1
           ;;
        -f|--file )
           file="$2"
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

#!/bin/bash

if [[ "$file" == "" ]]
then
    echo "usage: source validate_sgrna_list.sh [-f|--file] sgrnalist.txt"
    echo "       optional argument: --two-col-format (specify if input file is in two col format and requires conversion)"
    echo "file not specified. exiting."
    exit 1
fi

# quality check for metadata text file

# requirement

# 1. file exists

if [[ ! -e "$file" ]]
then
    echo "file $file does not exist. exiting."
    exit 1
fi

# 2. file has more than 0 bytes

if [[ ! -s "$file" ]]
then
    echo "file $file exists, but is empty. exiting."
    exit 1
fi

# 3. file is text file

file_type=$(file $file)
if [[ ! "$file_type" =~ ASCII ]]
then
    echo "file $file exists, but is not text (ASCII) file type. exiting."
    echo "file type: $file_type"
    exit 1
fi

# 4. file does not contain CRLF (windows newline)

if [[ "$(grep -U $'\015' $file 2>/dev/null)" != "" ]]
then
    echo "bad format for $file (CRLF detected -- convert to Unix format). exiting."
    exit 1
fi

# 5. file is rectangular tsv with header

ncol=$(head -n1 $file | tr '\t' '\n' | wc -l)
field_counts=($(gawk 'BEGIN{FS="\t"}{print NF}' $file))
flag=FALSE
for i in ${!field_counts[@]}
do
    field_count=${field_counts[$i]}
    if [[ "$field_count" != "$ncol" ]]
    then
        echo "row $i has $field_count fields ($ncol columns defined in the header)."
        flag=TRUE
    fi
done
if [[ "$flag" == "TRUE" ]]
then
    echo "exiting."
    exit 1
fi

# 6. file has exactly 3 col or exactly 2, depending on parameter
field_count=($(gawk 'BEGIN{FS="\t"}(FNR==1){print NF}' $file))

if [[ "$two_col_format" == "FALSE" ]]
then
    # typical
    if [[ $field_count != 3 ]]    
    then
        echo "file $file contains $field_count columns, where 3 are expected. exiting."
        exit 1
    fi
else
    # unusual
    if [[ $field_count != 2 ]]    
    then
        echo "file $file contains $field_count columns, where 1 are expected. exiting."
        exit 1
    fi
fi
