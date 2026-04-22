#!/bin/bash

#SBATCH --job-name="amr_NF_driver"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --partition=prod-compute,prod-compute-mem

# AMR Nextflow Pipeline - NCAH Modules Wrapper
# USAGE:
#   1. Navigate to your working directory where your data files are located.
#   2. Run this script with sbatch, providing --input_dir and --outdir flags.
#
# EXAMPLE:
#   sbatch /path/to/cluster_amr_nextflow_ncah_modules.sh \
#       --input_dir /path/to/my/data \
#       --outdir /path/to/results
#
#   # With custom MLST lookup file and Kraken database
#   sbatch /path/to/cluster_amr_nextflow_ncah_modules.sh \
#       --input_dir /path/to/my/data \
#       --outdir /path/to/results \
#       --mlst_lookup_json /path/to/custom/lookup.json \
#       --kraken_db /project/bioinformatic_databases/databases/kraken/k2_standard_20260226
#
# NOTES:
#   - This script can be run from any directory
#   - Pipeline location is fixed at: ${HOME}/git/gitlab/amr_nextflow/
#   - MLST lookup defaults to: ${HOME}/git/gitlab/amr_nextflow/bin/lookup_genome_size.json
#   - Uses environment modules via ncah_modules profile
#   - Kraken database defaults to: /project/bioinformatic_databases/databases/kraken/kraken2

set -euo pipefail

# Clean environment and load required modules
module purge
module load slurm
module load openjdk/17.0.8.1_1
module load nextflow/24.04.4

# Pipeline location (fixed)
PIPELINE_ROOT="${HOME}/git/gitlab/amr_nextflow"

# Default MLST lookup location
DEFAULT_MLST_LOOKUP="${PIPELINE_ROOT}/bin/lookup_genome_size.json"

# Verify pipeline exists
if [ ! -f "${PIPELINE_ROOT}/amr_nextflow.nf" ]; then
    echo "ERROR: amr_nextflow.nf not found at ${PIPELINE_ROOT}"
    echo "       Please ensure the pipeline is cloned to: ${HOME}/git/gitlab/amr_nextflow/"
    exit 1
fi

# Parse command line arguments
INPUT_DIR=""
OUTPUT_DIR=""
MLST_LOOKUP=""
NEXTFLOW_ARGS=()

while [ "$#" -gt 0 ]; do
    case "$1" in
        --input_dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --outdir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --mlst_lookup_json)
            MLST_LOOKUP="$2"
            shift 2
            ;;
        *)
            NEXTFLOW_ARGS+=("$1")
            shift
            ;;
    esac
done

# Validate required arguments
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: Both --input_dir and --outdir flags are mandatory."
    echo "Usage: sbatch $0 --input_dir /path/to/data --outdir /path/to/results [--mlst_lookup_json /path/to/lookup.json] [--kraken_db /path/to/db]"
    exit 1
fi

if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Set MLST lookup to default if not provided
if [ -z "$MLST_LOOKUP" ]; then
    MLST_LOOKUP="$DEFAULT_MLST_LOOKUP"
    echo "Using default MLST lookup: $MLST_LOOKUP"
fi

# Verify MLST lookup file exists
if [ ! -f "$MLST_LOOKUP" ]; then
    echo "WARNING: MLST lookup file not found at: $MLST_LOOKUP"
    echo "         Pipeline will use its configured default or may fail if required."
fi

# Execute pipeline
echo "Creating output directory and changing location..."
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
echo "Launch directory is now: $(pwd)"

echo ""
echo "Starting Nextflow run with environment modules..."
echo "Pipeline location: ${PIPELINE_ROOT}"
echo "Input directory: $(realpath "$INPUT_DIR")"
echo "Output directory: $(pwd)"
echo "MLST lookup: $MLST_LOOKUP"

nextflow run "${PIPELINE_ROOT}/amr_nextflow.nf" \
    -profile ncah_modules \
    --input_files "$(realpath "$INPUT_DIR")/*.{fastq.gz,fq.gz,fa,fasta}" \
    --mlst_lookup_json "$MLST_LOOKUP" \
    "${NEXTFLOW_ARGS[@]}"

echo ""
echo "Pipeline complete! Results are in $(pwd)"