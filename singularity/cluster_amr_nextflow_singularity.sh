#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1 # Number of CPUs to allow, 2 will use full resources on one node.
#SBATCH --cpus-per-task=48
#SBATCH --job-name="amr_driver"
#SBATCH --partition=cpu_compute
#SBATCH --export=NONE

# this file needs to be edited for each environment!
# edit sbatch params
# # edit module load calls
# # edit kraken_db path

# module load calls
module load openjdk/18.0.2.1
module load nextflow/22.10.7
module load singularity/3.8.5

#package files
for i in *.fast*; do
    n=`echo $i | sed 's/_.*//' | sed 's/\..*//'`
    echo "n is : $n"
    mkdir -p $n
    mv $i $n/
done

#get directory where files are
file_dir="$(pwd)/*/*"

# echo ${file_dir}
nextflow run ${HOME}/git/gitlab/amr_nextflow/amr_nextflow.nf -profile singularity_nbaf --files "${file_dir}" $@

# Remove empty work directory after Nextflow completes
# rm -rf work/