#!/bin/bash
 
#SBATCH --job-name=container_build
#SBATCH --output=./container_build__%j.out
#SBATCH --error=./container_build__%j.err
#SBATCH --account="aap mr scicomp hpc cnah users" 
 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --partition=scicomp-compute

module load apptainer/1.1.9

# Set up cache directory
export APPTAINER_CACHEDIR="${HOME}/.apptainer_cache"
mkdir -p $APPTAINER_CACHEDIR

# Clean up any old repository code to ensure a fresh clone
echo "Preparing to clone internal USDA repository..."
if [ -d "amr_nextflow" ]; then
    rm -rf amr_nextflow
fi

# Clone the repository directly to the host machine
GIT_SSH_COMMAND="ssh -o UserKnownHostsFile=/dev/null -o StrictHostKeyChecking=no" git clone git@arsiaam0vgit10.usda.net:Cameron.Norris/amr_nextflow.git

# Build the container
echo "Starting Apptainer build..."
apptainer build --fakeroot --force amr_nextflow.sif amr_nextflow.def

echo "Build complete."