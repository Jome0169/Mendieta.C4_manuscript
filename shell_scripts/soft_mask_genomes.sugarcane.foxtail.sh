#!/bin/bash
#SBATCH --job-name=soft_mask       # Job name
#SBATCH --partition=schmitz_p               # Partition (queue) name, i.e., highmem_p
#SBATCH --ntasks=1                          # Run a single task	
#SBATCH --cpus-per-task=24                  # Number of CPU cores per task
#SBATCH --mem=100gb                          # Job memory request
#SBATCH --time=96:00:00                     # Time limit hrs:min:sec
#SBATCH --output=soft_mask.%j.out         # Standard output log
#SBATCH --error=soft_mask.%j.err          # Standard error log
#SBATCH --mail-type=END,FAIL                # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=john.mendieta@uga.edu   # Where to send mail
#SBATCH --export=NONE                       # do not load any env variables to compute node

cd $SLURM_SUBMIT_DIR

source /apps/lmod/lmod/init/zsh

ml RepeatMasker/4.0.9-p2-foss-2019b-Perl-5.30.0

export OMP_NUM_THREADS=24

#RepeatMasker -pa 24 -lib maizeTE02052020.fa -e rmblast -q -noisy -xsmall -gff Zm-B73-REFERENCE-NAM-5.0.fa

#BuildDatabase -name foxtail_millet -engine ncbi Pmiliaceum.fa
#BuildDatabase -name sugarcane -engine ncbi Sspontaneum.fa

#RepeatModeler -database foxtail_millet -pa 6
#RepeatModeler -database sugarcane -pa 6


RepeatMasker -pa 24 -lib RM_16877.ThuNov42157362021/foxtail_repeat_library.fa -e rmblast -q -noisy -xsmall -gff Pmiliaceum.fa
RepeatMasker -pa 24 -lib RM_50429.SatNov61730042021/sugarcange_repeat_library.fa -e rmblast -q -noisy -xsmall -gff Sspontaneum.fa

