#!/bin/bash

#SBATCH --account=phd-apworley42
#SBATCH --cpus-per-task=32

module load gnu9
module load mpich

mpirun -n 1 ./gol.exe > out1.csv
mpirun -n 2 ./gol.exe > out2.csv
mpirun -n 4 ./gol.exe > out4.csv
mpirun -n 8 ./gol.exe > out8.csv
mpirun -n 16 ./gol.exe > out16.csv
mpirun -n 32 ./gol.exe > out32.csv
