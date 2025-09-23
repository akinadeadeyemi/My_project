#!/bin/bash
#SBATCH --job-name=skin_downsampled_correlation_N30_job
#SBATCH --output=skin_downsampled_correlation_N30_job.out
#SBATCH --error=skin_downsampled_correlation_N30_job.err
#SBATCH --mem=60gb
#SBATCH --time=15:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4


### Loading  R
module load r/4.4.0


### specifying the data path to be used for the correlation analysis
# final_blood_merged_betas="~/Phd_data/final_blood_merged_betas.csv"
# final_blood_samples="~/Phd_data/final_blood_samples.csv"



Rscript skin_downsampled_correlation_N30.R
