# Instructions for running via Singularity container

## Description
Use these instructions to set up Singularity container for AMR pipeline on HPC system. Currently working on NCAH and NBAF HPC systems


## Requirements
HPC needs to have Singularity and Nextflow installed. Development was done on a CentOS 7 system with Singularity 3.8.5 and Nextflow 22.10.6.

## Installation
Log into HPC 
1. Create a directory to store the container and Nextflow workflow
   `mkdir -p ${HOME}/git/gitlab/`
2. Clone the git repository (.ssh key must be set up)
`cd ${HOME}/git/gitlab/ && git clone git@arsiaam0vgit10.usda.net:Cameron.Norris/amr_nextflow.git`
3. Navigate to singularity directory and build container
   
   <u>If NCAH HPC:</u>

   `cd ${HOME}/git/gitlab/amr_nextflow/singularity/ && srun singularity build --fakeroot --tmpdir /local/$USER/ amr_nextflow.sif amr_nextflow.def`
   
   <u>If NBAF HPC:</u>

   `cd ${HOME}/git/gitlab/amr_nextflow/singularity/ && srun -p cpu_compute singularity build --fakeroot amr_nextflow.sif amr_nextflow.def`


## Usage
*** CURRENTLY CONFIGURED FOR NBAF HPC ***
1. Navigate to a directory containing any of these 3 file types:
   * paired-end .fastq.gz files
   * Single .fastq.gz file
   * .fasta assembly files
   
2. Call the `cluster_amr_nextflow_singularity.sh` script to organize into sample directories & execute Nextflow workflow.
To alter abricate coverage/depth parameters, you can pass `--abricate_depth` or `--abricate_coverage` parameters to this call. Defaults are 0 and 75 respectively.

   `sbatch ${HOME}/git/gitlab/amr_nextflow/singularity/cluster_amr_nextflow_singularity.sh`

   --or--

      `sbatch ${HOME}/git/gitlab/amr_nextflow/singularity/cluster_amr_nextflow_singularity.sh --abricate_depth 50 --abricate_coverage 90`

## Results
Results will collect in sample-specific folders. If you want to review workflow performance, each submission generates an html file that can be found in the original file directory. This report will also alert of any process failures.

## Development Notes
Need to implement config and sbatch script for running on local machine. Needs to query system cpu count and memory to set up nextflow.config file. 

## Support
Let me know if you have any questions or suggestions. Cameron.Norris@usda.gov


