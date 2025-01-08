#!/bin/bash

#SBATCH --job-name=shear-parallel         # create a short name for your job
#SBATCH --nodes=2                # node count
#SBATCH --ntasks=2               # total number of tasks across all nodes
#SBATCH --cpus-per-task=40        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=200G
#SBATCH --time=00:59:00          # total run time limit (HH:MM:SS)
#SBATCH -o /home/pa2937/logs/shear/shear.%A.out
#SBATCH -e /home/pa2937/logs/shear/shear.%A.err


########################################################################
#Can use SBATCH --mem-per-cpu=10G         # memory per cpu-core (4G is default)
#Can use SBATCH --mem=210G
# Can use SBATCH --mem-per-cpu=20G

module purge
module load anaconda3/2024.2
module load openmpi/gcc/4.1.1/64

export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=40
conda activate wlblends
cd /home/pa2937/2pt/wl-blends/2pt/

srun python mpi_pix_test.py

