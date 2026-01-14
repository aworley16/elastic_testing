#!/bin/bash

#SBATCH --account=phd-apworley42
#SBATCH --cpus-per-task=16

module load gnu9
module load mpich

mpirun -n 1 ./gol.exe
mpirun -n 2 ./gol.exe
mpirun -n 3 ./gol.exe
mpirun -n 4 ./gol.exe
mpirun -n 5 ./gol.exe
mpirun -n 6 ./gol.exe
mpirun -n 7 ./gol.exe
mpirun -n 8 ./gol.exe
mpirun -n 16 ./gol.exe