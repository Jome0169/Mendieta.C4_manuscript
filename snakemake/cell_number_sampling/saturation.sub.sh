#!/bin/bash
#SBATCH --job-name=FP_sat       # Job name
#SBATCH --partition=schmitz_hm_p               # Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task	
#SBATCH --cpus-per-task=20                  # Number of CPU cores per task
#SBATCH --mem=400gb                          # Job memory request
#SBATCH --time=96:00:00                     # Time limit hrs:min:sec
#SBATCH --output=fp_sat.%j.out         # Standard output log
#SBATCH --error=fp_sat.%j.err          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=john.mendieta@uga.edu   # Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node

cd $SLURM_SUBMIT_DIR

source /apps/lmod/lmod/init/zsh

ml BEDTools
ml SAMtools 


ml Anaconda3
source activate snakemake_6
snakemake -rps FP_saturation_analysis.cell_number.snake --use-conda --cores 20
