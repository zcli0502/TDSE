#!/bin/bash

#SBATCH --job-name=He4D_zcli
#SBATCH --partition=cpu
#SBATCH --mail-type=end
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH -n 4
#SBATCH --ntasks-per-node=1
#SBATCH --exclusive
##SBATCH --time=00:01:00 


OMP_NUM_THREADS=16
mpirun  ./he cfg_laser.cfg
