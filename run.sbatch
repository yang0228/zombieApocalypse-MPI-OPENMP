#!/bin/bash
#SBATCH --job-name="MPI-BS"
#SBATCH --nodes=2
#SBATCH --time=0-0:20
srun --ntasks-per-node=1 ./stage2