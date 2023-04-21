#!/bin/bash

# Shell script to submit an AF2 run consisting of two depending jobs.
# Run this script directly, don't invoke it via `sbatch`!

set -e

JOBID1=$(sbatch --parsable jobscript-alphafold2-step_1-msa.sh)
JOBID2=$(sbatch --parsable --dependency=afterok:${JOBID1} --deadline=now+2weeks jobscript-alphafold2-step_2-prediction.sh)

echo "Submitted jobs"
echo "    ${JOBID1} (MSA on CPU)"
echo "    ${JOBID2} (prediction on GPU)"
