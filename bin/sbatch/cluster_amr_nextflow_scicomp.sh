#!/bin/bash

#SBATCH --job-name="AMR_NF_driver"
#SBATCH --account="aap mr scicomp hpc cnah users" 

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --partition=scicomp-compute,scicomp-high-memory


module load apptainer/1.1.9

module purge
module load slurm
module load openjdk/17.0.11_9
module load nextflow/25.10.4
module load apptainer/1.1.9


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
nextflow run /project/training/pipeline_development/amr_nextflow/amr_nextflow/amr_nextflow.nf -profile scicomp_apptainer --files "${file_dir}" $@

# Remove empty work directory after Nextflow completes
rm -rf work/
