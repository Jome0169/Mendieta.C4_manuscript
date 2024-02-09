#!/bin/bash

##----------------------------------------------##
##          SLURM submission properties         ##
##----------------------------------------------##

#SBATCH --partition=highmem_p
#SBATCH --job-name=map_sciATAC_os_leaf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=7-00:00:00
#SBATCH --mem=120g
#SBATCH --output=logs_map_sciATAC_os_leaf.%j.log
#SBATCH --error=logs_map_sciATAC_os_leaf.%j.err

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
ref=00.data/Osativa_323_v7.0.fa
threads=24


declare -a rep_1_array=( "Rice_P1" )

bwa index ${ref}

for i in ${rep_1_array[@]}; do 
    echo "Working on file rep1" ${i} ; 

    # run BWA rep1
    bwa mem -M -t $threads $ref $fastq/${i}_R1.fastq.gz $fastq/${i}_R2.fastq.gz | samtools view -hbSq 10 | samtools sort - | perl util_scripts/modify_BC_flag.pl - | samtools view -bSh - > ${i}.mpq10.bam

done

#for i in ${rep_2_array[@]}; do 
#    echo "Working on file rep1" ${i} ; 
#
#    # run BWA rep1
#    bwa mem -M -t $threads $ref $fastq/rep2/${i}_R1.fastq.gz $fastq/rep2/${i}_R2.fastq.gz | samtools view -hbSq 30 | samtools sort - | perl modify_BC_flag.pl - | samtools view -bSh - > ${i}.rep2.mpq30.bam
#
#done


