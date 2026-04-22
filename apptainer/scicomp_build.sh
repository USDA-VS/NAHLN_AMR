#!/bin/bash
 
#SBATCH --job-name=container_build
#SBATCH --output=./container_build_%j.out
#SBATCH --error=./container_build_%j.err
#SBATCH --account="aap mr scicomp hpc cnah users" 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --partition=scicomp-compute

module load apptainer/1.1.9

export APPTAINER_CACHEDIR="${HOME}/.apptainer_cache"
export TMPDIR="/tmp/apptainer_${USER}_${SLURM_JOB_ID}"
export APPTAINER_TMPDIR="${TMPDIR}"

mkdir -p "${APPTAINER_CACHEDIR}"
mkdir -p "${TMPDIR}"

ls -ld "${APPTAINER_CACHEDIR}"
ls -ld "${TMPDIR}"

export APPTAINERENV_SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}

echo "Starting build with ${SLURM_CPUS_PER_TASK} CPUs..."
apptainer build --fakeroot --force amr_nextflow.sif amr_nextflow.def

echo "Build complete - woot woot."
ls -lh amr_nextflow.sif

rm -rf "${TMPDIR}"