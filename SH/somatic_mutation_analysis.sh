# Commands to analyse the identified somatic mutations

## Get transcription binding motifs in the somatic mutation regions and randomly selected regions
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/sm_analysis/smlist_20052022/
motifs=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/annotations/TF_binding_motifs/20220519132126_JASPAR2022_combined_matrices_20358_meme.txt

cd $CWD
bsub -q normal -n 1  -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo fimo_sm.log -eo fimo_sm.err "
    fimo -o fimo_sm $motifs sm_reg.fa
"
bsub -q normal -n 1  -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo fimo_sample.log -eo fimo_sample.err "
    fimo -o fimo_sample $motifs gen_samp.fa
"
