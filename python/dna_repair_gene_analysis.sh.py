"""
Functions to check whether the L1 vs L2 mutation rate differences arises because of
altered gene expression of DNA repair genes
"""
import pandas as pd
import numpy as np

# <editor-fold desc="Check DE differences in DNA repair genes">
# Get (cleaned) Thaliana gene ids
idmap = pd.read_table('TAIR2UniprotMapping.txt', header=None)
idmap[0] = [i.split('-')[0] for i in idmap[0]]
geneids = pd.read_table('DNA_repair_genes.txt', header=None)
geneids.drop_duplicates([0], inplace=True)
geneids.drop_duplicates([2], inplace=True)

# get rows with IDs
dbids = [s for s in geneids[2] if re.search(r'AT.G', str(s))]
geneids = geneids.loc[~geneids[2].isin(dbids)]

geneids[0] = geneids[0].apply(str.upper)
dbids2 = [s for s in geneids[0] if re.search(r'AT.G', str(s))]
geneids = geneids.loc[~geneids[0].isin(dbids2)]

dbids3 = [s for s in geneids[1] if re.search(r'AGI_LocusCode:AT.G', str(s))]
geneids = geneids.loc[~geneids[1].isin(dbids3)]

locusids = [s for s in geneids[1] if re.search(r'TAIR:locus:', str(s))]
geneids = geneids.loc[~geneids[1].isin(locusids)]

uniids = [s for s in geneids[1] if re.search(r'UniProtKB:', str(s))]
geneids = geneids.loc[~geneids[1].isin(uniids)]

assert geneids.shape[0] == 0

# Clean IDs
dbids = [d.split()[0] for d in dbids]
dbids2 = [d.split()[0] for d in dbids2]
dbids3 = [d.split(':')[1] for d in dbids3]
locusids = [d.split(':', maxsplit=1)[1] for d in locusids]
uniids = [d.split(':')[1] for d in uniids]

ldbid = list(idmap.loc[idmap[1].isin(locusids), 2].unique())
udbid = list(idmap.loc[idmap[0].isin(uniids), 2].unique())

alldbids = list(set(dbids + dbids2 + dbids3 + ldbid + udbid))

with open('TAIR_DNA_repair_gene_ids.txt', 'w') as fout:
    fout.write('\n'.join(alldbids))
# </editor-fold>

# <editor-fold desc='Get orthologues DNA repair genes in currot'>
# Done using python/find_athal_ortho_genes_in_cur.py:408
# # </editor-fold>

# <editor-fold desc='DNA repair genes overexpressed in Thaliana SAM'>
"""
Use list of single cell cluster marker genes from  https://www.sciencedirect.com/science/article/pii/S1534580721001611
to check if any of the DNA repair genes in TAIR_DNA_repair_gene_ids.txt are overexpressed in any of the SAM and leaf primordial clusters 
"""
cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/annotations/'
fin = '1-s2.0-S1534580721001611-mmc2.xlsx'
# From shell
# hometools xls2csv 1-s2.0-S1534580721001611-mmc2.xlsx 1-s2.0-S1534580721001611-mmc2_cluster-enriched_genes.txt -s 'cluster-enriched genes'
dnarep = [l.strip() for l in open(f'{cwd}/TAIR_DNA_repair_gene_ids.txt')]
clusters = {'Mesophyll': [0, 3, 4, 16],
            'shoot meristem': [1, 13],
            'epidermial': [2, 12, 14, 18],
            'proliferating': [5, 9, 17, 19],
            'vascular': [6, 7, 10, 22],
            'gaurd': [11],
            'companion': [15],
            'shoot endodermis': [20],
            'undefined': [8, 21]}
clsters2 = {v1: k for k, v in clusters.items() for v1 in v}

degenes = pd.read_csv(f'{cwd}/1-s2.0-S1534580721001611-mmc2_cluster-enriched_genes.txt')
degenes = degenes['avg_logFC p_val_adj cluster gene Description'.split()]
degenes = degenes.loc[degenes['gene'].isin(dnarep)]
degenes['type'] = [clsters2[i] for i in degenes['cluster']]
degenes.to_csv(f'{cwd}/dna_repair_DE_genes_in_thaliana_SAM_scRAN.txt', index=False, header=False, sep='\t')
# This identified 31 DNA repair genes that were present/distributed in different lineages. That is,
# there was no obvious enrichment of DNA repair genes in a layer (or cell type)
# </editor-fold>

# TODO: Analyse in R XXXXXXXXX

# TODO: Visualise



