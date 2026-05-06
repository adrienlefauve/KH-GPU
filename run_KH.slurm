#!/bin/bash
#SBATCH --job-name=KH
# Replace the account name with your own SLURM allocation
#SBATCH --account=YOUR_SLURM_ACCOUNT
# Replace the partition name with the GPU partition on your cluster
#SBATCH --partition=YOUR_GPU_PARTITION
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --time=16:00:00
#SBATCH --output=slurm_%j.log  # captures stdout/stderr as a backup


# Replace this with the absolute path to your local copy of the repo
cd /path/to/julia_simulation_code

# Show allocation for debugging
echo "Node: $(hostname)"
echo "Cores allocated: $SLURM_CPUS_PER_TASK"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"

# Set Julia threads to match CPUs allocated
export JULIA_NUM_THREADS=16

# Start the run
timestamp=$(date +%Y%m%d_%H%M%S)

julia start_KH_from_params.jl 2>&1 | tee "output_${timestamp}.txt"
