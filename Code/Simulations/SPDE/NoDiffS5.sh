#!/bin/bash

#$ -N RAs5ND
#$ -q math
#$ -pe openmp 64
#$ -cwd            		# run the job out of the current directory
#$ -m beas
#$ -ckpt blcr
#$ -o output/
#$ -e output/
module load MATLAB/r2015a
mkdir -p output
./driver 'No Diffusion S5' '1000' '.0001' '100' '5' '9' 'Flat' '0' '0' '0' '1' '1' '0'