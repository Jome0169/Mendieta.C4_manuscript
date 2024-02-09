#!/bin/bash

#SBATCH --partition=schmitz_hm_p
#SBATCH --job-name=process_dups_sciATAC
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
ml picard/2.21.6-Java-11
ml SAMtools 
# input files
#sci1=Sorghum_leaf.sciATAC_rep1.mq10.bam
#sci2=Sorghum_leaf.sciATAC_rep2.mq10.bam

declare -a rep_1_array=( ZP_P11_G029_10x ZP_P1_G029_10x ZP_P2_G029_10x ZP_P3_G029_10x )

# functions
doCall(){

	# input
	base=$1
	rep=$2
	threads=$3

	 remove duplicates
	echo "removing dups - $base ..."
	java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
		MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
		REMOVE_DUPLICATES=true \
		METRICS_FILE=$base.30mpq.metrics \
		I=$base.mpq30.bam \
		O=$base.mq30.BC.rmdup.bam \
		BARCODE_TAG=BC \
		ASSUME_SORT_ORDER=coordinate

    source activate perl 
	# fix barcodes
	echo "fixing barcodes and removing multi-mapped reads ..."
	perl ./util_scripts/fixBC.pl $base.mq30.BC.rmdup.bam $rep | samtools view -bhS - > $base.mq30.BC.rmdup.mm.bam

    source activate python3_dev
	# make Tn5 bed files
	echo "making Tn5 bed files ..."
	#perl ./util_scripts/makeTn5bed.pl $base.mq10.BC.rmdup.mm.bam | sort -k1,1 -k2,2n - > $base.tn5.bed
	samtools view $base.mq30.BC.rmdup.mm.bam | grep -v "XA:Z" | python makeTn5bed.py -sam - | sort -k1,1 -k2,2n - > $base.mpq30.tn5.bed

	# remove non-unique Tn5 sites
	uniq $base.mpq30.tn5.bed > $base.mpq30.unique.tn5.bed

}
export -f doCall

# run pipeline
for i in ${rep_1_array[@]}; do 
    echo "Working on file rep1" ${i} ; 

    doCall ${i} ${i} $threads
done

