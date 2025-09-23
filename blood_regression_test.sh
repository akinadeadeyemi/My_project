#!/bin/bash
#SBATCH --job-name=blood_regression_wt_lifespan_job
#SBATCH --output=blood_regression_wt_lifespan_job.out
#SBATCH --error=blood_regression_wt_lifespan_job.err
#SBATCH --mem=60gb
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4


### Loading  R
module load r/4.4.0


### specifying the data path to be used for the correlation analysis
# final_blood_merged_betas="~/Phd_data/final_blood_merged_betas.csv"
# final_blood_samples="~/Phd_data/final_blood_samples.csv"



Rscript blood_regression_test_wt_lifespan.R
