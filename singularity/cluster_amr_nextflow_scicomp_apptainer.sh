#!/bin/bash
 
#SBATCH --job-name="amr_NF_driver"
#SBATCH --account="aap mr scicomp hpc cnah users" 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --partition=scicomp-compute,scicomp-high-memory

# Clean environment and load required modules
module purge
module load slurm
module load openjdk/17.0.11_9
module load nextflow/25.10.4
module load apptainer/1.1.9

# Usage:
#   Assumes input is a collection of fastq files in the present working directory 
#   which is then made into a directory structure by this wrapper script.

echo "Packaging fastq files into directories..."
for i in *.fast*; do
    # Safety check: ensure the glob actually matched files
    if [ -f "$i" ]; then
        # Extract prefix (removes everything after the first underscore or dot)
        n=$(echo "$i" | sed 's/_.*//' | sed 's/\..*//')
        echo "n is : $n"
        mkdir -p "$n"
        mv "$i" "$n"/
    fi
done

# Get directory where files are now located
file_dir="$(pwd)/*/*"

# Run the pipeline
echo "Starting Nextflow run..."
nextflow run ${HOME}/git/gitlab/amr_nextflow/amr_nextflow.nf -profile scicomp_apptainer --files "${file_dir}" $@

# Remove empty work directory after Nextflow completes
# rm -rf work/