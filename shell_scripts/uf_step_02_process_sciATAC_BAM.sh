#!/bin/bash

#SBATCH --partition=schmitz_hm_p
#SBATCH --job-name=uf_process_dups_sciATAC
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
ml purge 
ml Anaconda3
ml picard/2.21.6-Java-11
ml SAMtools


# input files
#sci1=Sorghum_leaf.sciATAC_rep1.mq10.bam
#sci2=Sorghum_leaf.sciATAC_rep2.mq10.bam

#Input Array
declare -a Uro_sample_array=( Uro_P1 Uro_P2 Uro_P3 Uro_P4 Uro_P5 Uro_P6 Uro_P7 Uro_P8 Uro_T1 Uro_T2 )


# functions
doCall(){

	# input
	base=$1
	rep=$2
	threads=$3

	# remove duplicates
	#echo "removing dups - $base ..."
    #java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
    #    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    #    REMOVE_DUPLICATES=true \
    #    METRICS_FILE=$base.metrics \
    #    I=$base.mpq10.bam \
    #    O=$base.mq10.BC.rmdup.bam \
    #    BARCODE_TAG=BC \
    #    ASSUME_SORT_ORDER=coordinate

    source activate perl 
	# fix barcodes
	echo "fixing barcodes and removing multi-mapped reads ..."
	perl ./util_scripts/fixBC.pl $base.mq10.BC.rmdup.bam $rep | samtools view -bhS - > $base.mq10.BC.rmdup.mm.bam



	# make Tn5 bed files
	echo "making Tn5 bed files ..."
    source activate python3_dev
	samtools view $base.mq10.BC.rmdup.mm.bam | python makeTn5bed.py -sam - | sort -k1,1 -k2,2n - > $base.tn5.bed


	# remove non-unique Tn5 sites
	uniq $base.tn5.bed > $base.unique.tn5.bed

}
export -f doCall

# run pipeline
#doCall Sorghum_leaf.sciATAC_rep1 Sorghum_leaf.sciATAC_rep1 $threads
#doCall Sorghum_leaf.sciATAC_rep2 Sorghum_leaf.sciATAC_rep2 $threads



for i in ${Uro_sample_array[@]}; do 
    echo "Working on file rep1" ${i} ; 

    doCall ${i}.rep1 ${i}.rep1 $threads
    doCall ${i}.rep2 ${i}.rep2 $threads

done




