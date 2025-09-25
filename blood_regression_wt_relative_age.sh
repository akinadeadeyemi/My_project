#!/bin/bash
#SBATCH --job-name=blood_regression_relative_age_job
#SBATCH --output=blood_regression_relative_age_job.out
#SBATCH --error=blood_regression_relative_age_job.err
#SBATCH --mem=60gb
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4


### Loading  R
module load r/4.5.0

export R_LIBS_USER=/home/aakinad/R/x86_64-pc-linux-gnu-library/4.4

Rscript blood_regression_wt_relative_age.R
