#!/bin/bash
#SBATCH --job-name=skin_regression_wt_max_lifespan_job
#SBATCH --output=skin_regression_wt_max_lifespan_job.out
#SBATCH --error=skin_regression_wt_max_lifespan_job.err
#SBATCH --mem=60gb
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4


### Loading  R
module load r/4.4.0


Rscript skin_regression_max_lifespan_test.R
