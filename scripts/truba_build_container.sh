#!/bin/bash
# Build the Apptainer image that every later phase runs inside.
#
# TRUBA forbids installing conda/pip onto /arf (truba.md): those installs create
# hundreds of thousands of small files and degrade the shared Lustre filesystem.
# The sanctioned route is copy-open-write-close-use on a central base image.
#
#   nohup ./scripts/truba_build_container.sh > logs/container_build.log 2>&1 &
set -euo pipefail
ROOT=/arf/scratch/suysal/mtcovmap
BASE=/arf/sw/containers/miniconda3/miniconda3-container.sif
MM_URL=https://micro.mamba.pm/api/micromamba/linux-64/latest

cd "$ROOT"
# 1. copy
[ -f base.sif ] || cp "$BASE" base.sif
# 2. open
rm -rf sandbox
apptainer build --sandbox sandbox base.sif
# 3. write. conda install breaks under --fakeroot (signature-verification plugin
# conflict), so micromamba does the installing; the container may lack curl.
# A fresh prefix, not the base env: the base image ships a partially-populated
# conda package cache and installing into it dies with "Cannot find a valid
# extracted directory cache".
apptainer exec --writable --fakeroot sandbox bash -c '
  set -e
  mkdir -p /opt/bin
  wget -qO- '"$MM_URL"' | tar -xj -C /opt/bin --strip-components=1 bin/micromamba
  export MAMBA_ROOT_PREFIX=/opt/mamba
  /opt/bin/micromamba create -y -p /opt/mtcov -c conda-forge -c bioconda \
    python=3.11 pandas pyarrow scipy statsmodels matplotlib pyyaml pysam \
    samtools hisat2 star haplogrep
  /opt/bin/micromamba clean -y --all
  echo "export PATH=/opt/mtcov/bin:\$PATH" > /.singularity.d/env/99-mtcov.sh
'
# 4. close
rm -f mtcovmap.sif
apptainer build mtcovmap.sif sandbox
rm -rf sandbox base.sif
# 5. use - prove it before anything depends on it
apptainer exec --bind /arf mtcovmap.sif bash -c '
  python -c "import pandas, pysam, yaml, scipy, statsmodels, matplotlib; print(\"py ok\")"
  samtools --version | head -1
  STAR --version
  hisat2 --version | head -1
  haplogrep --version 2>&1 | head -1
'
echo "container ready: $ROOT/mtcovmap.sif"
