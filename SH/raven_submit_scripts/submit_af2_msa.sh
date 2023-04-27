sbatch -J falist.1 \
    -o output_%A_%a.txt \
    -e error_%A_%a.txt \
    /raven/u/mgoel/apricot/scripts/SH/raven_submit_scripts/jobscript-alphafold2-step_1-msa.sh \
    /raven/u/mgoel/apricot/cur_protein/fa.list.1.txt