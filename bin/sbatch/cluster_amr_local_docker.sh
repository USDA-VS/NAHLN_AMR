#!/bin/bash

# AMR Nextflow Pipeline - Local Docker Wrapper
# USAGE:
#   1. Navigate to the root directory of the cloned amr_nextflow pipeline.
#   2. Run this script, providing the --input_dir flag to point to your data.
#
# EXAMPLE:
#   ./bin/cluster_amr_local_docker.sh --input_dir /path/to/my/data

set -euo pipefail

# store root
PIPELINE_ROOT=$(pwd)

# make sure user is where they need to be
if [ ! -f "${PIPELINE_ROOT}/amr_nextflow.nf" ]; then
  echo "ERROR: amr_nextflow.nf not found in the current directory."
  echo "Please run this script from the root of the amr_nextflow pipeline repository."
  exit 1
fi

# check for tools, nextflow and docker
if ! command -v docker &> /dev/null; then
  echo "ERROR: Docker is not installed or not in your PATH."
  echo "       Please install Docker Desktop and ensure it is running."
  exit 1
fi
if ! command -v nextflow &> /dev/null; then
  echo "ERROR: Nextflow is not installed or not in your PATH."
  echo "       Please install Nextflow (see https://www.nextflow.io/docs/latest/getstarted.html)."
  exit 1
fi

# parse cmd line args
INPUT_DIR=""
OUTPUT_DIR=""
NEXTFLOW_ARGS=()

while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --outdir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        *)
            NEXTFLOW_ARGS+=("$1")
            shift
            ;;
    esac
done

# Validate Arguments
if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: Both --input_dir and --outdir flags are mandatory."
    echo "Usage: $0 --input_dir /path/to/data --outdir /path/to/results"
    exit 1
fi
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory not found: $INPUT_DIR"
    exit 1
fi

# Execute Pipeline
echo "Creating output directory and changing location..."
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
echo "Launch directory is now: $(pwd)"

echo ""
echo "Starting Nextflow run with Docker..."

nextflow run "${PIPELINE_ROOT}/amr_nextflow.nf" \
    -profile local_docker \
    --input_files "$(realpath "$INPUT_DIR")/*.{fastq.gz,fq.gz,fa,fasta}" \
    "${NEXTFLOW_ARGS[@]}"

echo ""
echo "Pipeline complete! Results are in $(pwd)"