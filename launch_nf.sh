#!/bin/bash
#SBATCH --partition cascade
#SBATCH --mem 5G
#SBATCH -c 1
#SBATCH -t 1:00:00

WORKFLOW=$1
CONFIG=$2

# Use a conda environment where you have installed Nextflow
# (may not be needed if you have installed it in a different way)
conda activate nextflow
module load Nextflow/23.04.2
module load Java/11.0.18
module load Anaconda3/2023.07-2
nextflow -C ${CONFIG} run ${WORKFLOW}  -with-trace -resume
