#!/bin/bash
#SBATCH --job-name=shear-subset         # create a short name for your job
#SBATCH --nodes=1                # node count
#SBATCH --ntasks=1               # total number of tasks across all nodes
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem-per-cpu=30G         # memory per cpu-core (4G is default)
#SBATCH --time=02:00:00          # total run time limit (HH:MM:SS)
#SBATCH -o /home/pa2937/logs/shear/shear.%A.out
#SBATCH -e /home/pa2937/logs/shear/shear.%A.err

module purge
module load anaconda3/2023.9
conda activate wlblends

cd /home/pa2937/2pt/wl-blends/2pt/

python subset_shear.py
