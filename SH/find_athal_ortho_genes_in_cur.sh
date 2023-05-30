CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/annotations
TAIR10GFF=/srv/biodata/dep_mercier/grp_schneeberger/data/Athal/TAIR10/TAIR10_GFF3_genes.gff
TAIR10CDS=/srv/biodata/dep_mercier/grp_schneeberger/data/Athal/TAIR10/TAIR10_cds_20101214
TAIR10PROT=/srv/biodata/dep_mercier/grp_schneeberger/data/Athal/TAIR10/TAIR10_pep_20101214
cd $CWD
cut -d',' -f 2 athal_cell_specific_genes.csv > athal_cell_specific_genes.geneids
# AT2G3830 is not a valid TAIR10 id and is therefore manually filtered out

grep -f athal_cell_specific_genes.geneids $TAIR10CDS \
| cut -d ' '  -f1 \
| sed 's/>//g' \
| sort > athal_cell_specific_genes.cdsids
hometools getchr -F athal_cell_specific_genes.cdsids -o athal_cell_specific_genes.cds.fa $TAIR10CDS
hometools getchr -F athal_cell_specific_genes.cdsids -o athal_cell_specific_genes.prot.fa $TAIR10PROT


# Blast against the CDS sequence of the CUR cds
curcds=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.pasa_out.CDS.fasta
blastn -query athal_cell_specific_genes.cds.fa \
     -task blastn \
     -db $curcds \
     -out cds.oblast \
     -outfmt 7 \
     -max_target_seqs 100 \
     -max_hsps 100 \
     -perc_identity 70 \
     -evalue 0.0001 \
     -num_threads 45 \
     >> blastall_blastn.log

# Blast against the prot sequence of the CUR prot
curprot=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.pasa_out.prot.fasta
blastp -query athal_cell_specific_genes.prot.fa \
     -task blastp \
     -db $curprot \
     -out prot.oblast \
     -outfmt 7 \
     -max_target_seqs 100 \
     -max_hsps 100 \
     -evalue 0.0001 \
     -num_threads 45 \
     >> blastall_blastn.log

# Run OrthoFinder
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/annotations/orthology
nohup orthofinder -t 40 -a 40 -f cur_athal &
# Get cur genes orthologous to A.thaliana marker genes
find_athal_ortho_genes_in_cur.py -> filter_orthologs()
#geneids="/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/annotations/athal_cell_specific_genes.geneids"
#orthogroups="/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/annotations/orthology/cur_athal/OrthoFinder/Results_Sep29/Orthogroups/Orthogroups.txt"
#python -c "
#import sys;
#sys.path.insert(0, '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/');
#from find_athal_ortho_genes_in_cur import filter_orthologs;
#filter_orthologs(\"$geneids\", \"$orthogroups\")
#" > cur_ortho_genes.txt
#
## Ortho Gene Functions
#grep -Ff cur_ortho_genes.txt /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/function/cur.pasa_out.prot.fasta.eggnog_mapper.tsv > cur_ortho_genes.func.txt

################################################################################
## Find cluster marker genes using thaliana leaf atlas
################################################################################
# In the paper (https://academic.oup.com/plcell/article/33/3/511/6067477?login=true)
# they list different marker genes for different cell types. I downloaded there
# list and then we can check the expression of "sequence orthologous" genes for
# the RP scRNA seq data

# Commands to find orthologous genes in find_athal_ortho_genes_in_cur.py -> get_cluster_marker_orthologs()


################################################################################
## Find orthologs between thaliana and apricot using structural similarity

### Create protein structure database
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/protein_structure/
cd $cwd
#### Athal structure DB
foldseek createdb athal_proteins/*pdb.gz athalprodb --threads 32
foldseek createindex athalprodb athalproindex --threads 32
#### P_armeniaca protein database
foldseek createdb prunus_armeniaca/*tar parmeniacaprodb --threads 32
foldseek createindex parmeniacaprodb parmeniacaproindex --threads 32

# TAIR UNIPROT ID MAPPING :  /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/protein_structure/TAIR2UniprotMapping-JAN2023.txt
# Example alignment
foldseek easy-search ~/Downloads/ranked_0.pdb /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/protein_structure/parmeniacaprodb out2.txt .