#!/bin/bash

#SBATCH --job-name=array
#SBATCH --output=logs/rr_array_%A_%a.out
#SBATCH --error=logs/rr_array_%A_%a.err
#SBATCH --array=1-3
#SBATCH --time=24:00:00
#SBATCH --partition=cdonnat
#SBATCH --mem=10G
#SBATCH --account=pi-cdonnat

# Print the task id.
echo "My SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
echo "My SLURM_ARRAY_JOB_ID: " $SLURM_ARRAY_JOB_ID
# Add lines here to run your computations
job_id=$SLURM_ARRAY_JOB_ID
module load R/4.2.0


result_file="new_normalized_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
echo "Result file is ${result_file}"
cd $SCRATCH/$USER/CCAR3/
Rscript experiments/simulations/simu_rr.R $SLURM_ARRAY_TASK_ID $result_file $1 $2 $3 $4 $5
#1: n
#2: strength theta
#3: p
#4: r
#5: q
