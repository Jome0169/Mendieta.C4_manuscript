#!/bin/bash
#SBATCH --job-name=RunHal       # Job name
#SBATCH --partition=batch               # Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task	
#SBATCH --cpus-per-task=5                  # Number of CPU cores per task
#SBATCH --mem=100gb                          # Job memory request
#SBATCH --time=120:00:00                     # Time limit hrs:min:sec
#SBATCH --output=ntr.%j.out         # Standard output log
#SBATCH --error=ntr.%j.err          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=john.mendieta@uga.edu   # Where to send mail

cd $SLURM_SUBMIT_DIR

source /apps/lmod/lmod/init/zsh

ml HAL

hal2maf evolverPlants_all_genomes.18_plants.hal evolverPlants_5_genomes.increased_gap.blocklen1000_maxgap50.V2.maf --refGenome Zm-B73  --noAncestors --noDupes --onlyOrthologs --maxBlockLen 1000 --maxRefGap 50 --targetGenomes Pmiliaceum,Sbicolor,Osativa,Ufusca
