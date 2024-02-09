input_name=${1}
tissue=${2}
base_name=${3}

source activate python3_dev

rg ${tissue} ${input_name} | awk 'BEGIN {OFS=FS="\t"}; {print $1,$2,$3,$4,$5,$6}' > ${base_name}.ortho.markers.bed
python process_markers.visualize.py -i ${base_name}.ortho.markers.bed  -o ${base_name}.ortho.visualize.bed -hd
python process_markers.annotate_v3.py -header -tis ${tissue} -marker ${base_name}.ortho.markers.bed -base ${base_name}.annotate
