#!/bin/bash

#SBATCH --partition=highmem_p
#SBATCH --job-name=pm_process_dups_sciATAC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=30
#SBATCH --time=7-00:00:00
#SBATCH --mem=490g
#SBATCH --output=logs_dedup_sciATAC.%j.log
#SBATCH --error=logs_dedup_sciATAC.%j.err

# set env
cd $SLURM_SUBMIT_DIR
source /apps/lmod/lmod/init/zsh


# threads
threads=30

# load modules
ml Anaconda3
ml SAMtools 
ml picard/2.21.6-Java-11

# input files
#sci1=Sorghum_leaf.sciATAC_rep1.mq10.bam
#sci2=Sorghum_leaf.sciATAC_rep2.mq10.bam

declare -a rep_1_array=( Proso_B1_G017_10x Proso_B1_G021_10x Proso_B9_G021_10x Proso_C1_G017_10x Proso_C9_G017_10x Proso_C9_G021_10x Proso_P1_G017_10x Proso_P2_G021_10x Proso_P3_G021_10x Proso_P8_G017_10x )

declare -a rep_2_array=( Proso_R2_P10_G021_10x Proso_R2_P11_G021_10x Proso_R2_P12_G021_10x Proso_R2_P13_G021_10x Proso_R2_P1_G017_10x Proso_R2_P1_G021_10x Proso_R2_P20_G017_10x Proso_R2_P8_G021_10x Proso_R2_P9_G021_10x )




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
		I=$base.mpq30.bam \
		O=$base.mq30.BC.rmdup.bam \
		BARCODE_TAG=BC \
		ASSUME_SORT_ORDER=coordinate

    source activate perl 
	# fix barcodes
	echo "fixing barcodes and removing multi-mapped reads ..."
	perl ./util_scripts/fixBC.pl $base.mq30.BC.rmdup.bam $rep | samtools view -bhS - > $base.mq30.BC.rmdup.mm.bam

	# make Tn5 bed files
    source activate python3_dev
	samtools view $base.mq30.BC.rmdup.mm.bam | python makeTn5bed.py -sam - | sort -k1,1 -k2,2n - > $base.mpq30.tn5.bed

	# remove non-unique Tn5 sites
	#uniq $base.tn5.bed > $base.unique.tn5.bed

	uniq $base.mpq30.tn5.bed > $base.unique.mpq30.tn5.bed

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



