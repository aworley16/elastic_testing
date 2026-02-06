#!/bin/bash

#SBATCH --account=phd-apworley42
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --exclusive

module load gnu9
module load mpich

mpirun -n 1 ./gol.exe 1024 6000  >> singlee.csv
mpirun -n 2 ./gol.exe 1024 6000  >> singlee.csv
mpirun -n 4 ./gol.exe 1024 6000  >> singlee.csv
mpirun -n 8 ./gol.exe 1024 6000  >> singlee.csv
mpirun -n 16 ./gol.exe 1024 6000 >> singlee.csv

#mixed phases
mpirun -n 1 ./gol.exe 1024 2000 4out124.csv 3 1 2 >> 1_2_4e.csv

mpirun -n 1 ./gol.exe 1024 3000 out18.csv 2 1 8 >> 18e.csv

mpirun -n 4 ./gol.exe 1024 3000 out48.csv 2 4 8 >> 48.csv

mpirun -n 1 ./gol.exe 1024 2000 out124.csv 3 1 2 4>> 1_2_4p.csv
