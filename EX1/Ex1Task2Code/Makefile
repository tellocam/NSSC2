# requirements on ubuntu
# sudo apt-get build-essentials
# sudo apt-get install openmpi-bin openmpi-common libopenmpi-dev

# required modules on cluster
# module load mpi/openmpi-x86_64
# module load pmi/pmix-x86_64

all: 
	mpic++ stenciljacobiMPI1D.cpp -o jacobiMPI -O3 -pedantic -Wall -std=c++14 -march=native -ffast-math -lpthread -DUSEMPI