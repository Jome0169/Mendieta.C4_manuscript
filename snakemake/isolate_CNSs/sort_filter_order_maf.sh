#!/usr/bin/env bash
set -euxo pipefail
original_maf_file=${1}
base_file_name=${original_maf_file%.*}


sed s/\.Zm_b73//g ${original_maf_file} > ${base_file_name}.fixed_chrom.maf

mafOrder ${base_file_name}.fixed_chrom.maf order_list.txt ${base_file_name}.ordered.maf

maf-sort ${base_file_name}.ordered.maf > ${base_file_name}.ordered.sorted.maf

mafFilter -speciesFilter=Zm.files.txt ${base_file_name}.ordered.sorted.maf > ${base_file_name}.ordered.sorted.filtered.maf

maf_parse -g Zm-B73.scored_20bp_windows.final.CNS.sorted.no_exons.modified.bed ${base_file_name}.ordered.sorted.filtered.maf > ${base_file_name}.intersecting_CNS.maf

maffilter input.file=${base_file_name}.intersecting_CNS.maf input.file.compression=none "maf.filter=Merge(species=(Zm-B73,Pmiliaceum,Sbicolor,Osativa,Ufusca), dist_max=20), OutputAlignments(file=${base_file_name}.intersecting_CNS.blocks.clustal, mask=no, coordiantes=yes)"
