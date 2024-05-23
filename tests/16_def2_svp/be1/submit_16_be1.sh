#!/bin/bash

#SBATCH -t 1:00:00
#SBATCH -o 16_def2-svp.out
#SBATCH -e 16_def2-svp.err
#SBATCH -c 8
#SBATCH -N 1
#SBATCH --partition=short

source /home/minsik/.bashrc
eval "$(conda shell.bash hook)"
conda activate split

ulimit -s unlimited
export OMP_STACKSIZE=500m
export OMP_NUM_THREADS=8

python3 ../../../parallel_split_code.py ./16.xyz be1 def2-svp 8 0 > output.log
