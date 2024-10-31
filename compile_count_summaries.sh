mode=fastq
ncol=8
first_sample=$(gawk '(FNR==2){print $2}' sample_metadata.txt)
cat \
	<(gawk \
		-v ncol=$ncol \
		'(FNR==1){printcol=ncol==""?NF:ncol; printf $1; for(i=2;i<=printcol;i++){printf OFS $i} printf ORS}' \
		output/${first_sample}_$mode/${first_sample}_$mode.countsummary.txt) \
	<(
		for samp in $(gawk '(FNR>1&&$0!~/^#/){print $2}' sample_metadata.txt)
		do
			gawk \
				-v ncol=$ncol \
				'BEGIN{FS=OFS="\t"}(FNR>1){printcol=ncol==""?NF:ncol; printf $1; for(i=2;i<=printcol;i++){printf OFS $i} printf ORS}' \
				output/${samp}_$mode/${samp}_$mode.countsummary.txt \
				2>/dev/null
		done
	)
