#!/bin/bash
#SBATCH --job-name=skin_phylo_regression_job
#SBATCH --output=skin_phylo_regression_job.out
#SBATCH --error=skin_phylo_regression_job.err
#SBATCH --mem=60gb
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4



### Correcting for inabilities to load phytools properly
# Loadin interactive shell environment
source ~/.bashrc  # or ~/.bash_profile if you're using macOS or zsh



### Loading  R
module load r/4.4.0


### specifying the data path to be used for the correlation analysis
# final_blood_merged_betas="~/Phd_data/final_blood_merged_betas.csv"
# final_blood_samples="~/Phd_data/final_blood_samples.csv"



Rscript skin_phylo_regression.R
