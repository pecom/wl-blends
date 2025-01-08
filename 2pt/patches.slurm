#!/bin/bash

#SBATCH --job-name=buzzard-patches         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=100G
#SBATCH --time=4:00:00          # total run time limit (HH:MM:SS)
#SBATCH -o /home/pa2937/logs/shear/patches.%A.out
#SBATCH -e /home/pa2937/logs/shear/patches.%A.err


########################################################################
#Can use #SBATCH --mem-per-cpu=10G         # memory per cpu-core (4G is default)

module purge
module load anaconda3/2024.2

# cd /home/pa2937/
conda activate wlblends
cd /home/pa2937/2pt/wl-blends/2pt/

python make_bigfile.py
# srun python mpi_shear.py
# srun python mpi_pix_test.py

