#!/bin/bash

#SBATCH --partition=schmitz_hm_p
#SBATCH --job-name=zm_process_dups_sciATAC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=3-00:00:00
#SBATCH --mem=490g
#SBATCH --output=logs_dedup_sciATAC.%j.log
#SBATCH --error=logs_dedup_sciATAC.%j.err

# set env
cd $SLURM_SUBMIT_DIR
source /apps/lmod/lmod/init/zsh


# threads
threads=10

# load modules
ml Anaconda3
ml SAMtools 
ml picard/2.21.6-Java-11

# input files
#sci1=Sorghum_leaf.sciATAC_rep1.mq10.bam
#sci2=Sorghum_leaf.sciATAC_rep2.mq10.bam


declare -a rep_1_array=( Zm_rep1_P10_10x Zm_rep1_P1_10x Zm_rep1_P2_10x Zm_rep1_P3_10x Zm_rep1_P4_10x Zm_rep1_P5_10x Zm_rep1_P6_10x Zm_rep1_P7_10x Zm_rep1_P8_10x Zm_rep1_P9_10x Mo17B73_P2_G031_10x )

#
declare -a rep_2_array=( Zm_rep2_P1_10x Zm_rep2_P2_10x Zm_rep2_P3_10x Zm_rep2_P4_10x Zm_rep2_P5_10x Zm_rep2_P6_10x Mo17B73_P2_G031_10x )



# functions
doCall(){

	# input
	base=$1
	rep=$2
	threads=$3

	## remove duplicates
	echo "removing dups - $base ..."
	java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
		REMOVE_DUPLICATES=true \
		METRICS_FILE=$base.metrics \
		I=$base.mpq10.bam \
		O=$base.mq10.BC.rmdup.bam \
		BARCODE_TAG=BC \
		ASSUME_SORT_ORDER=coordinate

    source activate perl 
	# fix barcodes
	echo "fixing barcodes and removing multi-mapped reads ..."
	perl ./util_scripts/fixBC.pl $base.mq10.BC.rmdup.bam $rep | samtools view -bhS - > $base.mq10.BC.rmdup.mm.bam

	# make Tn5 bed files
    source activate python3_dev
	samtools view $base.mq10.BC.rmdup.mm.bam | python util_scripts/makeTn5bed.py -sam - | sort -k1,1 -k2,2n - > $base.mpq10.tn5.bed

	# remove non-unique Tn5 sites
	#uniq $base.tn5.bed > $base.unique.tn5.bed

	uniq $base.mpq10.tn5.bed > $base.unique.mpq10.tn5.bed

}
export -f doCall

# run pipeline
#doCall Sorghum_leaf.sciATAC_rep2 Sorghum_leaf.sciATAC_rep2 $threads


for i in ${rep_1_array[@]}; do 
    echo "Working on file rep1" ${i} ; 

    doCall ${i}.rep1 ${i}.rep1 $threads
done

for i in ${rep_2_array[@]}; do 
    echo "Working on file rep2" ${i} ; 

    doCall ${i}.rep2 ${i}.rep2 $threads

done



