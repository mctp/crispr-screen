first_sample=$(gawk '(FNR==2){print $2}' sample_metadata.txt)
cat \
	<(awk '(FNR==1){print $0}' output/$first_sample/$first_sample.countsummary.txt) \
	<(
		for samp in $(gawk '(FNR>1&&$0!~/^#/){print $2}' sample_metadata.txt)
		do
			tail -n +2 output/$samp/$samp.countsummary.txt 2>/dev/null
		done
	)
