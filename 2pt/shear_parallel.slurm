#!/bin/bash

#SBATCH --job-name=shear-parallel         # create a short name for your job
#SBATCH --nodes=2                # node count
#SBATCH --ntasks=60               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1:00:00          # total run time limit (HH:MM:SS)
#SBATCH -o /home/pa2937/logs/shear/shear.%A.out
#SBATCH -e /home/pa2937/logs/shear/shear.%A.err


########################################################################
#Can use SBATCH --mem-per-cpu=10G         # memory per cpu-core (4G is default)
#Can use SBATCH --mem=210G

module purge
module load anaconda3/2024.2
module load openmpi/gcc/4.1.1/64

conda activate wlblends
cd /home/pa2937/2pt/wl-blends/2pt/

srun python mpi_pix_test.py

