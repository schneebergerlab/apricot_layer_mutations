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


# <editor-fold desc="Filtering genes close to SMs.">


# Initially developed for Jose /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data_sharing/functional_data_for_Jose_13_04_2023/filtering_scripts.sh
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/
anno=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/EVM_PASA/pasa_on_mancur/cur.pasa_out.sort.protein_coding.3utr.gff3
strmap=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/structurome/structural_orthologs_cur_athal.geneids.tsv
cd $cwd

window_size=1000
for s in wt_1 wt_19 mut_11_1 mut_15; do
    # Get mutation BED files
    grep -P "\t$s\t" all_sm_in_all_samples.manually_selected.cleaned.csv \
     | awk '{print $1"\t"$2-1"\t"$2}' \
     | sort | uniq > ${s}.sms.bed
    # Select genes from GFF within the window
    bedtools window -w $window_size -a ${s}.sms.bed -b $anno \
    | grep -P '\tgene\t' \
    | cut -f12 \
    | cut -d';' -f1 \
    | cut -d'=' -f2 \
    | sort | uniq  > ${s}.genes_with_sm_within_1kb.txt
    # Get Athal IDs for these genes
    grep -Ff ${s}.genes_with_sm_within_1kb.txt $strmap > ${s}.genes_with_sm_within_1kb_with_athal.txt
done


## Change gene ids to mRNA Ids
#sed 's/Gene/mRNA/g' selected_genes.txt > selected_mRNAs.txt
#
## Get functional for the selected mRNA. Can also pipe with "cut -F " to select specific columns from  cur.pasa_out.prot.fasta.eggnog_mapper.tsv
#grep -Ff selected_mRNAs.txt cur.pasa_out.prot.fasta.eggnog_mapper.tsv > selected_mRNA_function.txt
#
## Select flowering time related genes
#bedtools window -w $window_size -a mut.bed -b cur_athal_ortho_gene_annotations.bed | cut -f7 | sort | uniq  > selected_flowering_genes.txt

# </editor-fold>

# <editor-fold desc="Cluster WT and SM protein structures>
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/structurome/smeffect/
cd $cwd
# Create and index DB
foldseek createdb pdbs/ pdbdb --threads 32
foldseek createindex pdbdb tmp --threads 32
# Align and cluster
foldseek search pdbdb pdbdb aln tmp --alignment-type 1 -c 0.5 -a -e 0.001
foldseek clust pdbdb aln clu --cluster-mode 0
foldseek createtsv pdbdb pdbdb clu clu.tsv
foldseek convertalis --threads 40 pdbdb pdbdb aln search.m8

# </editor-fold>