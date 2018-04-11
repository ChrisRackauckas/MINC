#!/bin/bash

#$ -N RAs3ND
#$ -q math
#$ -pe openmp 64
#$ -cwd            		# run the job out of the current directory
#$ -m beas
#$ -ckpt blcr
#$ -o output/
#$ -e output/
module load MATLAB/r2015a
mkdir -p output
./driver 'No Diffusion S3' '1000' '.0001' '100' '3' '9' 'Flat' '0' '0' '0' '1' '1' '0'