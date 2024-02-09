cicero_files=${1}
marker_genes=${2}
output_prefix=${3}

bedtools intersect -a ${cicero_files} -b ${marker_genes} -wa -wb > ${output_prefix}.cicero_intersecting_markers.bed











