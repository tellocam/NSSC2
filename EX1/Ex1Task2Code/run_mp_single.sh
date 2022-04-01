#!/bin/bash
#SBATCH -J jacobiMPI
#SBATCH -N 1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#SBATCH --exclusive 
#SBATCH --time=00:10:00
if command -v sinfo  2>/dev/null # if on cluster
then
    module load mpi/openmpi-x86_64
    module load pmi/pmix-x86_64
    mpiproc=20
else  # if on local machine
    mpiproc=8
fi

mpirun -n $mpiproc ./jacobiMPI 128 1000
