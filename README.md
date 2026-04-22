# amr_nextflow

## Description
Antimicrobial resistance (AMR) pipeline using Nextflow with dependencies managed by Apptainer/Docker.

----

## General information

The AMR pipeline provides AMR identification, taxonomic classification, and genome assembly for next-generation sequening data and assemblies.

<p align="center">
  <img src="https://raw.githubusercontent.com/USDA-VS/NAHLN_AMR/main/docs/img/AMR.png" width="500" height="700">
</p>

----

## Dependencies

- [Nextflow >= v25.10.4](https://docs.seqera.io/nextflow/install)
- [Apptainer](https://apptainer.org/docs/admin/main/installation.html) (Required for HPC execution) OR [Docker Desktop](https://www.docker.com/products/docker-desktop/) (Required for Mac Laptop execution)
   - If using Apptainer, pull the Apptainer image (instructions below)
- [Kraken2 database](https://benlangmead.github.io/aws-indexes/k2) 

*Nextflow can be installed either [system-wide](https://docs.seqera.io/nextflow/install) or with [Conda](https://anaconda.org/channels/bioconda/packages/nextflow/overview).

------

## Run the AMR pipeline

### Quick Start

**The following `nextflow run...` commands must be ran from the GitHub repository location (i.e., /amr_nextflow/).**

**w/Docker on Laptop**
```bash
nextflow run amr_nextflow.nf \
-profile local_docker \
--input_files "/path/to/input/files/*.{fastq.gz,fq.gz,fa,fasta}"\
 --kraken_db /path/to/kraken2/database
```

**w/Apptainer on HPC w/o SLURM**
```bash
nextflow run amr_nextflow.nf \
    -profile hpc_apptainer \
    --input_files "/home/user.name/test_files/*.{fastq.gz,fq.gz,fa,fasta}" \
    --outdir "/home/user.name/results-hpc_apptainer" \
    --kraken_db /path/to/databases/kraken/k2_standard_20260226
```

If SLURM is require on a server, you can wrap the above in an sbatch script. This will submit all job under a single SLURM job ID.

Alternatively, you may generate a new profile for your specific HPC, if you are not using the NCAH HPC or the NBAF SciComp HPC, following these profiles as examples.

SLURM sbatch wrapper script example:
```bash
#!/bin/bash

#SBATCH --job-name="amr_solo"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --export=NONE
#SBATCH --partition=prod-compute,prod-compute-mem

set -euo pipefail

# Clean environment and load required modules
module purge
module load slurm
module load openjdk/17.0.8.1_1
module load nextflow/25.10.4
module load apptainer/1.1.9

nextflow run amr_nextflow.nf \
    -profile hpc_apptainer \
    --input_files "/home/user.name/test_files/*.{fastq.gz,fq.gz,fa,fasta}" \
    --outdir "/home/user.name/results-hpc_apptainer" \
    --kraken_db /path/to/databases/kraken/k2_standard_20260226
```
*Note the SLURM options, flags and partition names will differ on each HPC. Additionally, module versions may slightly differ and can be replaced by Conda or system-wide installations.*

------

## Command-line options

```bash
usage: nextflow run amr_nextflow.nf [-profile] [--input_files] [-outdir] [--kraken_db] [-resume]

General:
   --input_files        STR      Path to input directory containing input files (fastq.gz,fq.gz,fa,fasta)
   --outdir             STR      Path to output directory. (default: ./results)

   --abricate_depth     FLOAT    Minimum DNA %coverage of the resistance gene for a hit to be reported.
                                 This controls how much of the gene must be present in your assembly.
                                 (Corresponds to Abricate's --mincov flag).
                                 (default: 0)

   --abricate_identity  FLOAT    Minimum DNA %identity of the resistance gene for a hit to be reported.
                                 This controls how similar your assembly must be to the gene.
                                 (Corresponds to Abricate's --minid flag).
                                 (default: 75)

Databases:
   --kraken_db          STR      Path to preformatted Kraken2 database. MUST be downloaded. (default location: /kraken/)

Container:
   --amr_nf_container   PATH     Path to SIF apptainer image. Only required if manually supplying Apptainer image. (default: Automatic Docker download and conversion to Apptainer SIF)

Reporting:
   --logo               STR      Path to PNG logo for use in PDF reports. (deault: /doc/img/default.png)

Nextflow Options: 
    -resume                 Pipeline will resume from previous run if terminated
    -with-report            A single document which includes many useful metrics about a workflow execution
    -with-trace             Creates an execution tracing file that contains some useful information about each process.executed in your pipeline script
    -with-timeline          Render an HTML timeline for all processes executed in your pipeline
    -with-dag               creates a file containing a textual representation of the pipeline execution graph in the DOT format     
```

------


## Detailed Use-Cases for Mac Laptop, NCAH HPC and the NBAF SciComp HPC

## Nextflow Profiles for use on NCAH HPC, NBAF SciComp and a Mac laptop:

There are 3 separate Nextflow profiles:

1. `ncah_apptainer` ran by executing the `cluster_amr_nextflow_ncah_apptainer.sh` wrapper.
2. `scicomp_apptainer` ran by executing the `cluster_amr_nextflow_scicomp_apptainer.sh` wrapper.
3. `local_docker` ran by executing the `cluster_amr_nextflow_local_docker.sh` wrapper.

## Kraken databases

A [Kraken database](https://benlangmead.github.io/aws-indexes/k2) must be provided for taxonomic classification with Kraken2.

For example, on the NCAH HPC we have these:

 ```bash
ls /project/bioinformatic_databases/databases/kraken/

k2_standard_20260226
kraken2_std8
kraken2 #aka kraken2_host
kraken2_max
#etc
#etc
 ```

Download Kraken databases [here](https://benlangmead.github.io/aws-indexes/k2)

-----

## How to run on NCAH:

**These instrtuctions enable the use of the `cluster_amr_nextflow_ncah_apptainer.sh` wrapper script, specifically designed for the NCAH HPC.**

1. Make the directory structure in your ${HOME}, if it does not exist:

```bash
cd ${HOME}

mkdir -p git/gitlab/

cd git/gitlab/
```

2. Pull the GitLab repository and switch to the `upgrade/apptainer` branch:

```bash
git clone https://github.com/USDA-VS/NAHLN_AMR.git
```


3. Navigate to the data to be analyzed

*While the pipeline can be executed from (anywhere), we will execute in the location of the raw data being analyzed.*

Example:

Example data:

```bash
cd testing_dir

ls 

sample_1_R1.fastq.gz
sample_1_R2.fastq.gz
```

4. Run the following command, ensuring the input and output paths are specified correctly and that the Kraken database desired is used

```bash
sbatch ${HOME}/git/gitlab/amr_nextflow/bin/sbatch/cluster_amr_nextflow_ncah_apptainer.sh \
--input_dir /home/User.Name/testing_dir \
--outdir /home/User.Name/testing_dir/results \
--kraken_db /project/bioinformatic_databases/databases/kraken/k2_standard_20260226
```

*Kraken: The default is the large, core nt database, here we are using the smaller k2 standard database.*

*Paths: Full paths are used here for illustration and clarity. The pipeline does not require absolute paths.*

-----

## How to run on a Mac:

1. Download Docker Desktop:

**For USDA employees, you may need IT help, or just install it to your ${HOME}**

**Also, for USDA emplyees, you MUST run pipeline WITHIN the Docker Desktop App.**
**To open a teminal, click the ">_ terminal" button in the bottom right hand corner of the Docker Desktop App.**

[Download here](https://www.docker.com/products/docker-desktop/)

3. Install Conda/Anaconda and build a Conda environmetn with [Nextflow v25.10.4](https://anaconda.org/channels/bioconda/packages/nextflow/overview)

Instructions to install conda [here](https://docs.conda.io/projects/conda/en/stable/user-guide/install/index.html)

```bash
conda create -n nextflow-25 nextflow=25.10.4 -c conda-forge -y

conda activate nextflow-25
```
*[Nextflow >= v25.10.4](https://docs.seqera.io/nextflow/install) can also be install system wide*

3. Download the Kraken database of your choosing and note the location you download it to. 

Example where to download from the NCAH HPC:

```bash
/project/bioinformatic_databases/databases/kraken/kraken2_std8
```

*Kraken note: Kraken loads the database into memory. Thus, don't use a database bigger than your RAM.*
*The kraken2_std8 is a good, small database.*

Example download location:

```bash
/Users/user.name/databases/kraken2_std8
```

2. Create a directory to store the repo in ( does NOT have to be git/gitlab/ ), pull the repo, and switch to the `upgrade/simplify` branch:

Example:

```bash
cd  ${HOME}

git clone https://github.com/USDA-VS/NAHLN_AMR.git
```

3. Execute the pipeline FROM the git repo:

Example data:

```bash
ls /path/to/sample_data/

sample_1_R1.fastq.gz
sample_1_R2.fastq.gz
```

```bash
cd amr_pipeline

bash bin/sbatch/cluster_amr_local_docker.sh \
--input_dir /path/to/sample_data \
--outdir /path/to/sample_data/results \
--kraken_db /Users/user.name/databases/kraken2_std8
```

*Execution location: The local version of this pipeline MUST to be executed from the repo location.*

-----

## How to run on SciComp:

**These instrtuctions enable the use of the `cluster_amr_nextflow_scicomp_apptainer.sh` wrapper script, specifically designed for the NBAF SciComp HPC.**

1. Make the directory structure in your ${HOME}, if it does not exist:

```bash
cd ${HOME}

mkdir -p git/gitlab/

cd git/gitlab/
```

2. Pull the GitLab repository and switch to the `upgrade/apptainer` branch:

```bash
git clone https://github.com/USDA-VS/NAHLN_AMR.git
```

3. Navigate to the data to be analyzed

*While the pipeline can be executed from (anywhere), we will execute in the location of the raw data being analyzed.*

Example data:

```bash
cd testing_dir

ls 

sample_1_R1.fastq.gz
sample_1_R2.fastq.gz
```

4. Run the following command, ensuring the input and output paths are specified correctly and that the Kraken database desired is used

```bash
sbatch ${HOME}/git/gitlab/amr_nextflow/bin/sbatch/cluster_amr_nextflow_scicomp_apptainer.sh \
--input_dir /home/User.Name/testing_dir \
--outdir /home/User.Name/testing_dir/results \
--kraken_db /project/training/temp_databases/k2_standard_20260226

*Kraken: The default is the large, core nt database, here we are using the smaller k2 standard database.*

*Paths: Full paths are used here for illustration and clarity. The pipeline does not require absolute paths.*
