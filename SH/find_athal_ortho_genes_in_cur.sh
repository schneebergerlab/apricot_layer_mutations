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

#### Apricot protein database
## Select highest scoring PDB model for each mRNA
## python.find_athal_ortho_genes_in_cur.select_pdb_from_af2

foldseek createdb currot_pdbs/selected_pdb.tar curdb --threads 32
foldseek createindex curdb tmp --threads 32

# align currot to athal
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/structurome/
indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/protein_structure/
cd $cwd
bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=20000]" -M 25000 -oo athal_cur.log -eo athal_cur.err "
    foldseek search --threads 40 ${indir}/curdb ${indir}/athalprodb athal_cur_search_db tmp -c 0.5 -a -e 0.001
    foldseek convertalis --threads 40 ${indir}/curdb ${indir}/athalprodb athal_cur_search_db athal_cur_search.tsv --format-output \"query,target,fident,nident,alnlen,qcov,tcov,qstart,qend,tstart,tend,evalue,prob,bits,mismatch,gapopen\" --format-mode 4
    foldseek rbh --threads 40 ${indir}/curdb ${indir}/athalprodb athal_cur_rbh_db tmp -c 0.5 -a -e 0.001
    foldseek convertalis --threads 40 ${indir}/curdb ${indir}/athalprodb athal_cur_rbh_db athal_cur_rbh.tsv --format-output \"query,target,fident,nident,alnlen,qcov,tcov,qstart,qend,tstart,tend,evalue,prob,bits,mismatch,gapopen\" --format-mode 4
"


# Test script/commands
cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/tests/foldseek_test
foldseek createdb athal/ athalDB
foldseek createdb prunus/ prunusDB
foldseek createindex prunusDB tmp
foldseek createindex athalDB tmp
foldseek search prunusDB athalDB alnDB tmp -c 0.5 -a -e 0.001
foldseek convertalis prunusDB athalDB alnDB aln.tsv --format-output "query,target,fident,nident,alnlen,qcov,tcov,qstart,qend,tstart,tend,evalue,prob,bits,mismatch,gapopen" --format-mode 4
foldseek rbh prunusDB athalDB rbhDB tmp -c 0.5 -a -e 0.001
foldseek convertalis prunusDB athalDB rbhDB rbh.tsv --format-output "query,target,fident,nident,alnlen,qcov,tcov,qstart,qend,tstart,tend,evalue,prob,bits,mismatch,gapopen" --format-mode 4



foldseek result2msa prunusDB athalDB alnDB msa --msa-format-mode 6  # Generate MSA
foldseek unpackdb msa msa_output --unpack-suffix a3m --unpack-name-mode 0