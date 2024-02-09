#!/bin/bash

# vars
name=$1
peaks=$( ls *summits.bed )

# iterate
for i in ${peaks[@]}; do

	echo "$i"

	id=$( echo "$i" | cut -d'.' -f1-2 )	

	perl peak_calling_scripts/normalize_score.pl $i \
	| perl -ne 'chomp;my@col=split("\t",$_);
		$col[1] = $col[1] - 250;
		$col[2] = $col[2] + 250;
		if($col[1] >= 0){
			print "$col[0]\t$col[1]\t$col[2]\t$col[3]\t$col[0]_$col[1]_$col[2]_$col[3]\n";}' - > $id.temp
done

# merge peaks
cat *.temp | sort -k1,1 -k2,2n - > $1.temp2

# actually merge
bedtools merge -i $1.temp2 -c 5 -o collapse | perl peak_calling_scripts/selectNonOverlapping.pl - > $1.unique500bpPeaks.bed
rm *.temp *.temp2
