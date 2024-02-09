#!/bin/bash

##----------------------------------------------##
##          SLURM submission properties         ##
##----------------------------------------------##

#SBATCH --partition=highmem_p
#SBATCH --job-name=map_sciATAC_Uf_leaf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=7-00:00:00
#SBATCH --mem=120g
#SBATCH --output=logs_map_sciATAC_Uf_leaf.%j.log
#SBATCH --error=logs_map_sciATAC_Sb_leaf.%j.err

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
ref=00.data/Ufusca_669_v1.0.fa
threads=24

declare -a Uro_sample_array=( Uro_P1 Uro_P2 Uro_P3 Uro_P4 Uro_P5 Uro_P6 Uro_P7 Uro_P8 Uro_T1 Uro_T2 )


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


#bwa index ${ref}

for i in ${Uro_sample_array[@]}; do 
    echo "Working on file rep1" ${i} ; 

    # run BWA rep1
    bwa mem -M -t $threads $ref $fastq/rep1/${i}_10x_R1.fastq.gz $fastq/rep1/${i}_10x_R2.fastq.gz | samtools view -hbSq 10 | samtools sort - | perl modify_BC_flag.pl - | samtools view -bSh - > ${i}.rep1.mpq10.bam
 
    echo "Working on file rep2" ${i} ; 
    # run BWA rep1
    bwa mem -M -t $threads $ref $fastq/rep2/${i}_10x_R1.fastq.gz $fastq/rep2/${i}_10x_R2.fastq.gz | samtools view -hbSq 10 | samtools sort - | perl modify_BC_flag.pl - | samtools view -bSh - > ${i}.rep2.mpq10.bam
    
done


