#!/bin/bash

#$ -N RAcpu
#$ -q math
#$ -pe openmp 64
#$ -cwd            		# run the job out of the current directory
#$ -m beas
#$ -o output/
#$ -e output/
#$ -t 3-11
#$ -ckpt blcr
module load MATLAB/r2014b
mkdir -p output
./driver 'HPC Large Run, rndVar=5' 1e-8 5e6 500 5 $SGE_TASK_ID 1e-1 1 1