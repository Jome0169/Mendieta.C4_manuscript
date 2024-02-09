#!/bin/bash
#SBATCH --job-name=Demultiplex       # Job name
#SBATCH --partition=highmem_p               # Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task
#SBATCH --cpus-per-task=10                  # Number of CPU cores per task
#SBATCH --mem=400gb                          # Job memory request
#SBATCH --time=1:00:00                     # Time limit hrs:min:sec
#SBATCH --output=log.demulti.%j.out         # Standard output log
#SBATCH --error=log.demulti.%j.err          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=john.mendieta@uga.edu   # Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node

# set env
cd $SLURM_SUBMIT_DIR
source /apps/lmod/lmod/init/zsh

# load modules
ml SAMtools/1.10-GCC-8.3.0
ml Anaconda3

# bam
bamdir=/scratch/jpm73279/comparative_single_cell/02.QC_clustering/zea_mays/doublet_estimation

# vcf
vcfdir=/scratch/jpm73279/comparative_single_cell/02.QC_clustering/zea_mays/doublet_estimation
vcffile=$vcfdir/B73.Mo17.maize.biallelic.maf05.dp_filt.PHASED.vcf

# functions
runDemux(){

	# inputs
	bam=$1
	vcf=$2
	out=$3
	id=$4

	# sample info
	samples=/scratch/jpm73279/comparative_single_cell/02.QC_clustering/zea_mays/doublet_estimation/genotypes.txt
	meta=Zm_Mo17B73_P2_G031_10x.rep1_bc_counts.txt

	# mkdir	
	if [ ! -d $out ]; then
		mkdir $out
	fi

	# subset genotypes
	if [ ! -f $out/$id.INPUT.genotypes.txt ]; then

		if [ $id == "pool21" ] || [ $id == "pool22" ]; then
			awk -F'\t' -v poolID='pool1' '$2==poolID' $samples | cut -f1 - | uniq - > $out/$id.INPUT.genotypes.txt
		else
			awk -F'\t' -v poolID=$id '$2==poolID' $samples | cut -f1 - | uniq - > $out/$id.INPUT.genotypes.txt
		fi
	fi

	# subset barcodes
	if [ ! -f $out/$id.BARCODES.txt ]; then
		awk -F'\t' -v poolID=$id '$2==poolID' $meta | cut -f1 - | uniq - > $out/$id.BARCODES.txt
		sed -i 's/CB:/BC:/' $out/$id.BARCODES.txt
	fi
	
    conda activate popscle
	# run demuxlet -- change call rate to 0.9
	popscle demuxlet --sam $bam \
		--vcf $vcf \
		--out $out/$id.demuxlet \
		--sm-list $out/$id.INPUT.genotypes.txt \
		--group-list $out/$id.BARCODES.txt \
		--field GT \
		--min-MQ 10 \
		--min-callrate 1

}
export -f runDemux


runDemux Zm_Mo17B73_P2_G031_10x.rep1.mq10.BC.rmdup.mm.renamed_BC_CB.bam $vcffile B73_Mo16_doublet B73_Mo16_doublet

#pool1=$( seq 1 40 )
#for i in $pool1; do
#	runDemux $bamdir/pool$i.BC.mq10.rmdup.bam $vcffile $PWD/pool$i pool$i
#done
   
