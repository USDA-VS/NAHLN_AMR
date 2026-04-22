#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=2 # Number of CPUs to allow, 2 will use full resources on one node.
#SBATCH --job-name="amr_driver"
#SBATCH --partition=prod-compute,prod-compute-mem
#SBATCH --export=NONE

module purge
module load slurm
module load openjdk/17.0.8.1_1
module load nextflow/23.10.0

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
nextflow run ${HOME}/git/gitlab/amr_nextflow/amr_nextflow.nf -profile slurm_ncah --files "${file_dir}" $@

# Remove empty work directory after Nextflow completes
rm -rf work/