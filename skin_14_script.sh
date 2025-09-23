#!/bin/bash
#SBATCH --job-name=skin_sample_14_correlation_job
#SBATCH --output=skin_sample_14_correlation_job.out
#SBATCH --error=skin_sample_14_correlation_job.err
#SBATCH --mem=60gb
#SBATCH --time=15:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4


### Loading  R
module load r/4.4.0


### specifying the data path to be used for the correlation analysis
# final_blood_merged_betas="~/Phd_data/final_blood_merged_betas.csv"
# final_blood_samples="~/Phd_data/final_blood_samples.csv"



Rscript skin_sample_14_correlation.R




echo " Finished working on the correlation of all species"

echo "Thank you!"

