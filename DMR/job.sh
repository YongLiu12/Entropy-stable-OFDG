#!/bin/sh
#SBATCH -N 8
srun -n 512 ./all
