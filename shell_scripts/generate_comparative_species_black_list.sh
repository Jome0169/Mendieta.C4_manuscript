#!/bin/bash
#SBATCH --job-name=Collision_rate       # Job name
#SBATCH --partition=batch               # Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task	
#SBATCH --cpus-per-task=24                  # Number of CPU cores per task
#SBATCH --mem=100gb                          # Job memory request
#SBATCH --time=48:00:00                     # Time limit hrs:min:sec
#SBATCH --output=log_collision.%j.out         # Standard output log
#SBATCH --error=log_collision.%j.err          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=john.mendieta@uga.edu   # Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node

cd $SLURM_SUBMIT_DIR

source /apps/lmod/lmod/init/zsh

ml SAMtools 
ml BWA
ml BEDTools 
ml bioawk
threads=15
set -euxo pipefail

genome_1=Zm-B73-REFERENCE-NAM-5.0.fa
genome_2=Pmiliaceum.fa
hybrid_assembly="Zm_Proso.genome_combination.fa"

genome_1_base="Zm"
genome_2_base="Pm"

#bwa index ${genome_1}
#bwa index ${genome_2}
#bwa index ${hybrid_assembly}
#
#wgsim -N 50000000 -1 107 -2 107 -r 0 -R 0 -X 0 -e 0 -S 12345 ${genome_1} ${genome_1_base}.R1.fq ${genome_1_base}.R2.fq
#wgsim -N 50000000 -1 107 -2 107 -r 0 -R 0 -X 0 -e 0 -S 12345 ${genome_2} ${genome_2_base}.R1.fq ${genome_2_base}.R2.fq


#bioawk -c fastx '{print "@Zea_mays_"$name" "$comment"\n"$seq"\n+\n"$qual}' ${genome_1_base}.R1.fq > ${genome_1_base}.R1.renamed.fq 
#bioawk -c fastx '{print "@Zea_mays_"$name" "$comment"\n"$seq"\n+\n"$qual}' ${genome_1_base}.R2.fq > ${genome_1_base}.R2.renamed.fq
#bioawk -c fastx '{print "@Proso_millet_"$name" "$comment"\n"$seq"\n+\n"$qual}' ${genome_2_base}.R1.fq > ${genome_2_base}.R1.renamed.fq 
#bioawk -c fastx '{print "@Proso_millet_"$name" "$comment"\n"$seq"\n+\n"$qual}' ${genome_2_base}.R2.fq > ${genome_2_base}.R2.renamed.fq



#cat ${genome_1_base}.R1.renamed.fq ${genome_2_base}.R1.renamed.fq > merged_reads.${genome_1_base}.${genome_2_base}.R1.fq
#cat ${genome_1_base}.R2.renamed.fq ${genome_2_base}.R2.renamed.fq > merged_reads.${genome_1_base}.${genome_2_base}.R2.fq

bwa mem -M -t $threads $hybrid_assembly merged_reads.${genome_1_base}.${genome_2_base}.R1.fq merged_reads.${genome_1_base}.${genome_2_base}.R2.fq | samtools view -hbSq 30 | samtools sort - | samtools view -Sbh - > merged_genome.aligned_reads_from.${genome_1_base}.${genome_2_base}.mpq10.bam

bedtools bamtobed -i merged_genome.aligned_reads_from.${genome_1_base}.${genome_2_base}.mpq10.bam > merged_genome.aligned_reads_from.${genome_1_base}.${genome_2_base}.mpq10.bed


#Older Version Cross mapping genome 1 vs Genome 2 etc...
#bwa mem -M -t $threads $genome_1 ${genome_2_base}.R1.fq ${genome_2_base}.R2.fq | samtools view -hbSq 30 | samtools sort - | samtools view -Sbh - > ${genome_1_base}.aligned_reads_from.${genome_2_base}.mpq10.bam
#bwa mem -M -t $threads $genome_2 ${genome_1_base}.R1.fq ${genome_1_base}.R2.fq | samtools view -hbSq 30 | samtools sort - | samtools view -Sbh - > ${genome_2_base}.aligned_reads_from.${genome_1_base}.mpq10.bam
#
#
#bedtools bamtobed -i ${genome_1_base}.aligned_reads_from.${genome_2_base}.mpq10.bam > ${genome_1_base}.aligned_reads_from.${genome_2_base}.mpq10.bed
#bedtools bamtobed -i ${genome_2_base}.aligned_reads_from.${genome_1_base}.mpq10.bam > ${genome_2_base}.aligned_reads_from.${genome_1_base}.mpq10.bed




