#!/bin/bash
#SBATCH --job-name=blood_phylolm_model_job
#SBATCH --output=blood_phylolm_model_job.out
#SBATCH --error=blood_phylolm_model_job.err
#SBATCH --mem=60gb
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4


### Loading  R
module load r/4.4.0


Rscript blood_phylolm_model.R
