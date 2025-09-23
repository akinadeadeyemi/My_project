#!/bin/bash
#SBATCH --job-name=skin_phylolm_model_job
#SBATCH --output=skin_phylolm_model_job.out
#SBATCH --error=skin_phylolm_model_job.err
#SBATCH --mem=60gb
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4


### Loading  R
module load r/4.4.0


Rscript skin_phylolm_model.R
