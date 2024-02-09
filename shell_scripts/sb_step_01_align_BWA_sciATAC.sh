#!/bin/bash

##----------------------------------------------##
##          SLURM submission properties         ##
##----------------------------------------------##

#SBATCH --partition=highmem_p
#SBATCH --job-name=map_sciATAC_Sb_leaf
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=7-00:00:00
#SBATCH --mem=120g
#SBATCH --output=logs_map_sciATAC_Sb_root.%j.log
#SBATCH --error=logs_map_sciATAC_Sb_root.%j.err

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
ref=00.data/Sorghum_bicolor_var_BTx623.mainGenome.MtPt.fasta
threads=24

# rep1
A1=Sorghum_rep1_R1.fastq
A2=Sorghum_rep1_R2.fastq
Ai=Sorghum_rep1_index.fastq
Ai1=Sorghum_rep1_R1.BC.fastq
Ai2=Sorghum_rep1_R2.BC.fastq

# rep2
B1=Sorghum_rep2_R1.fastq
B2=Sorghum_rep2_R2.fastq
Bi=Sorghum_rep2_index.fastq
Bi1=Sorghum_rep2_R1.BC.fastq
Bi2=Sorghum_rep2_R2.BC.fastq

# attach barcodes
attachBC(){
	
	# read file
	read1=$1
	
	# index file
	readi=$2

	# output file
	out=$3

	# run
	umi_tools extract --bc-pattern=NNNNNNNNNNNNNNNNNNNNNNNNNNN --stdin=$readi --read2-in=$read1 --stdout=$out --read2-stdout

}
export -f

# attach barcodes
#attachBC $fastq/$A1 $fastq/$Ai $fastq/$Ai1
#attachBC $fastq/$A2 $fastq/$Ai $fastq/$Ai2
#attachBC $fastq/$B1 $fastq/$Bi $fastq/$Bi1
#attachBC $fastq/$B2 $fastq/$Bi $fastq/$Bi2


# run BWA
bwa mem -M -t $threads $ref $fastq/$Ai1 $fastq/$Ai2 | samtools view -hbSq 10 | samtools sort - | perl modify_BC_flag.pl - | samtools view -bSh - > Sorghum_leaf.sciATAC_rep1.mq10.bam

bwa mem -M -t $threads $ref $fastq/$Bi1 $fastq/$Bi2 | samtools view -hbSq 10 | samtools sort - | perl modify_BC_flag.pl - | samtools view -bSh - > Sorghum_leaf.sciATAC_rep2.mq10.bam
