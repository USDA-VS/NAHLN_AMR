#!/bin/bash
 
#SBATCH --job-name=container_build
 
#SBATCH --output=./container_build__%j.out
#SBATCH --output=./container_build__%j.err
 
#SBATCH --nodes=1
#SBATCH --partition=cpu_compute
#SBATCH --ntasks-per-node=20
 
module load singularity/3.8.5
singularity build --fakeroot --force amr_nextflow.sif amr_nextflow.def