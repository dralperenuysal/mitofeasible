#!/bin/bash
# Package the prepared sandbox into a .sif, on a compute node.
#
# squashfs compression of a ~3 GB sandbox exceeds the login node's CPU-time
# ulimit and dies with "signal: CPU time limit exceeded (core dumped)" after
# the packages are already installed. Only this last step needs a real node.
#
#   cd /arf/scratch/suysal/mtcovmap && sbatch <repo>/scripts/truba_package_container.slurm
#SBATCH --job-name=mtcov_sif
#SBATCH --partition=barbun
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=01:00:00
#SBATCH --output=/arf/scratch/suysal/mtcovmap/sif_%j.out
#SBATCH --error=/arf/scratch/suysal/mtcovmap/sif_%j.err

set -euo pipefail
ROOT=${MTCOV_ROOT:-/arf/scratch/suysal/mtcovmap}   # scratch: big intermediates
cd "$ROOT"
rm -f mtcovmap.sif
apptainer build mtcovmap.sif sandbox
rm -rf sandbox base.sif
# Prove the image before anything depends on it: an under-populated .sif fails
# silently at submit time and loudly hours later.
apptainer exec --bind /arf mtcovmap.sif bash -c '
  python -c "import pandas, pysam, yaml, scipy, statsmodels, matplotlib; print(\"py ok\")"
  samtools --version | head -1
  STAR --version
  hisat2 --version | head -1
  haplogrep --version 2>&1 | head -1
'
echo "container ready: $ROOT/mtcovmap.sif"
