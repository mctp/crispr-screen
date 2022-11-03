#!/bin/bash

# supply a fastq and read ID (and optionally a search sequence to highlight)
# ID is e.g. @A00839:131:HCNHJDRX2:1:2147:14145:16720

# usage: bash fsatq_record_by_id.sh <fastq.gz> <read_id> <optional: search_sequence>

FQ="$1"
RID="$2"
SRCH_SEQ="$3"

function revcomp() {
	if [[ -p /dev/stdin ]]
	then
		while read -r data
		do
			echo "$data" \
			| gawk \
				-v txtcase="$TXTCASE" \
                		-v nuctype="$NUCTYPE" \
		                'BEGIN{
		                        if(nuctype=="rna"){
                		                comp["u"]="a";comp["a"]="u"
                                		comp["U"]="A";comp["A"]="U"
		                        }else{
                		                comp["a"]="t";comp["t"]="a";
                                		comp["A"]="T";comp["T"]="A";
		                        }
                		        comp["c"]="g";comp["g"]="c"
		                        comp["C"]="G";comp["G"]="C"

                		        comp["r"]="y";comp["y"]="r"
		                        comp["R"]="Y";comp["Y"]="R"
		                }{
                		        seq=""
		                        for(i=length($0);i>0;i--){
                		                nuc=substr($0,i,1);
                                		if(comp[nuc]!=""){
		                                        nuc=comp[nuc]
                		                }
                                		seq = seq nuc
		                        }
                		        if(txtcase=="lower"){
                                		seq = tolower(seq)
		                        }else if(txtcase=="upper"){
                		                seq = toupper(seq)
		                        }
                		        print seq
		                }'
		done
	else
		echo "$1" \
		| gawk \
			-v txtcase="$TXTCASE" \
               		-v nuctype="$NUCTYPE" \
	                'BEGIN{
	                        if(nuctype=="rna"){
               		                comp["u"]="a";comp["a"]="u"
                               		comp["U"]="A";comp["A"]="U"
	                        }else{
               		                comp["a"]="t";comp["t"]="a";
                               		comp["A"]="T";comp["T"]="A";
	                        }
               		        comp["c"]="g";comp["g"]="c"
	                        comp["C"]="G";comp["G"]="C"
                	        comp["r"]="y";comp["y"]="r"
	                        comp["R"]="Y";comp["Y"]="R"
	                }{
               		        seq=""
	                        for(i=length($0);i>0;i--){
               		                nuc=substr($0,i,1);
                               		if(comp[nuc]!=""){
	                                        nuc=comp[nuc]
               		                }
                               		seq = seq nuc
	                        }
               		        if(txtcase=="lower"){
                               		seq = tolower(seq)
	                        }else if(txtcase=="upper"){
               		                seq = toupper(seq)
	                        }
               		        print seq
	                }'
	fi
}

if [[ "$(file $FQ)" =~ "ASCII text" ]]
then
	CAT="cat"
else
	CAT="zcat"
fi

if [[ "$SRCH_SEQ" == "" ]]
then
	$CAT $FQ \
	| gawk -v rid="$RID" 'BEGIN{flag=0; flag_ct= 0; rid= rid " "}{ if(flag){ if(flag_ct < 4){ print $0; flag_ct++}else{exit}} if($0~rid){ print $0; flag=1; flag_ct=1}}'
else
	$CAT $FQ \
	| gawk -v rid="$RID" 'BEGIN{flag=0; flag_ct= 0; rid= rid " "}{ if(flag){ if(flag_ct < 4){ print $0; flag_ct++}else{exit}} if($0~rid){ print $0; flag=1; flag_ct=1}}' \
	| grep -iE --color=always "$(echo $SRCH_SEQ | revcomp)|$SRCH_SEQ|$"

fi
