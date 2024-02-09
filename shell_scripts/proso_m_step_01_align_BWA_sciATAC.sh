#!/bin/bash

##----------------------------------------------##
##          SLURM submission properties         ##
##----------------------------------------------##

#SBATCH --partition=highmem_p
#SBATCH --job-name=map_sciATAC_Proso_leaf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=7-00:00:00
#SBATCH --mem=120g
#SBATCH --output=logs_map_sciATAC_Prosm_leaf.%j.log
#SBATCH --error=logs_map_sciATAC_Prosom_leaf.%j.err

# set env
cd $SLURM_SUBMIT_DIR
source /apps/lmod/lmod/init/zsh

# modules
ml BWA/0.7.17-GCC-8.3.0
ml UMI-tools/1.0.1-foss-2019b-Python-3.7.4
ml SAMtools 

# common variables
fastq=00.data/fastq
bamout=BAM
ref=00.data/Pmiliaceum.fa
threads=24


declare -a rep_1_array=( Proso_B1_G017_10x Proso_B1_G021_10x Proso_B9_G021_10x Proso_C1_G017_10x Proso_C9_G017_10x Proso_C9_G021_10x Proso_P1_G017_10x Proso_P2_G021_10x Proso_P3_G021_10x Proso_P8_G017_10x )

declare -a rep_2_array=( Proso_R2_P10_G021_10x Proso_R2_P11_G021_10x Proso_R2_P12_G021_10x Proso_R2_P13_G021_10x Proso_R2_P1_G017_10x Proso_R2_P1_G021_10x Proso_R2_P20_G017_10x Proso_R2_P8_G021_10x Proso_R2_P9_G021_10x )


# attach barcodes
#attachBC(){
#	
#	# read file
#	read1=$1
#	
#	# index file
#	readi=$2
#
#	# output file
#	out=$3
#
#	# run
#	umi_tools extract --bc-pattern=NNNNNNNNNNNNNNNNNNNNNNNNNNN --stdin=$readi --read2-in=$read1 --stdout=$out --read2-stdout
#
#}
#export -f
#
## attach barcodes
##attachBC $fastq/$A1 $fastq/$Ai $fastq/$Ai1
##attachBC $fastq/$A2 $fastq/$Ai $fastq/$Ai2
##attachBC $fastq/$B1 $fastq/$Bi $fastq/$Bi1
##attachBC $fastq/$B2 $fastq/$Bi $fastq/$Bi2
#


bwa index ${ref}

for i in ${rep_1_array[@]}; do 
    echo "Working on file rep1" ${i} ; 

    # run BWA rep1
    bwa mem -M -t $threads $ref $fastq/rep1/${i}_R1.fastq.gz $fastq/rep1/${i}_R2.fastq.gz | samtools view -hbSq 30 | samtools sort - | perl modify_BC_flag.pl - | samtools view -bSh - > ${i}.rep1.mpq30.bam

done

for i in ${rep_2_array[@]}; do 
    echo "Working on file rep1" ${i} ; 

    # run BWA rep1
    bwa mem -M -t $threads $ref $fastq/rep2/${i}_R1.fastq.gz $fastq/rep2/${i}_R2.fastq.gz | samtools view -hbSq 30 | samtools sort - | perl modify_BC_flag.pl - | samtools view -bSh - > ${i}.rep2.mpq30.bam

done


