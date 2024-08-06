#!/bin/bash
#SBATCH --partition=radiology,fn_long
#SBATCH --job-name=gre_simulation
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --time=01-00:00:00
#SBATCH --mem=256GB
#SBATCH --array=1-5

echo $HOSTNAME
echo $SLURM_ARRAY_TASK_ID
START="$(date +%s)"

/gpfs/scratch/asslaj01/julia-1.10.0/bin/julia --threads=40  \
../main_slurm.jl \
--config_file_path "example_config.yaml" \
--num_jobs 5 \
--array_num $SLURM_ARRAY_TASK_ID \
--filesave "../results/test/kout_$SLURM_ARRAY_TASK_ID" \

DURATION=$[ $(date +%s) - ${START} ]
echo "time taken for simulation[s]:"
echo ${DURATION}