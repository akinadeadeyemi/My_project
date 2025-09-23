#!/bin/bash
#SBATCH --job-name=skin_regression_relative_age_residuals_ID_job
#SBATCH --output=skin_regression_relative_age_residuals_ID_job.out
#SBATCH --error=skin_regression_relative_age_residuals_ID_job.err
#SBATCH --mem=60gb
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4


### Loading  R
module load r/4.5.0

export R_LIBS_USER=/home/aakinad/R/x86_64-pc-linux-gnu-library/4.4

Rscript skin_regression_relative_age_wt_residuals_ID_test.R
