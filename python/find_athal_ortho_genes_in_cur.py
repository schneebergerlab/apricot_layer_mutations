import sys

import numpy as np

sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/software/hometools')

def filter_orthologs():
    """
        Get Currot genes that are orthologous to the cell-type specific marker genes from AThal.
        geneids: file contatining list of genes whose orthologs are to be selected
        orthogroups: path to Orthogroups.txt file generated by OrthoFinder
    """

    # <editor-fold desc="Defins imports">
    from collections import deque, defaultdict
    import pandas as pd
    import re
    from copy import deepcopy
    # </editor-fold>


    # <editor-fold desc="Define constants and defaults">
    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/annotations/'
    orthogroups = "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/annotations/orthology/cur_athal/OrthoFinder/Results_Sep29/Orthogroups/Orthogroups.txt"
    kkgenefin = "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/annotations/athal_cell_specific_genes.csv" # Cell-type gene ids identified by Kristin
    lagenefin = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/annotations/leaf_cell_type_markers/tpc.00751.2020-s03.csv'     # Leaf atlas cell type markers
    la2genefin = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/annotations/leaf_cell_type_markers2/table1.csv'     # Leaf atlas cell type markers 2
    la2layerfin = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/annotations/leaf_cell_type_markers2/Supplemental_dataset_3.xlsx'
    la2epidermfin = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/annotations/leaf_cell_type_markers2/Supplemental_dataset_6.xlsx'
    uniprot2tairfin = '/srv/biodata/dep_mercier/grp_schneeberger/data/Athal/TAIR10.1/TAIR2UniprotMapping-JAN2023.txt'
    curmrnafin = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/protein_structure/currot_pdbs/all_pdbs.txt'
    allathaluni = '/srv/biodata/dep_mercier/grp_schneeberger/data/Athal/TAIR10.1/uniprot-compressed_true_download_true_fields_accession_2Creviewed_2C-2023.06.07-13.09.58.00.tsv'
    strorthofin = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/structurome/structural_orthologs_cur_athal.geneids.tsv'
    # </editor-fold>


    # <editor-fold desc="Get orthologs for KK marker genes">
    genes = pd.read_table(kkgenefin, header=None, sep=',')
    genes.columns = 'cell_type athal_gene'.split()
    genesids = set(genes.athal_gene)
    ortho = defaultdict(deque)
    for l in open(orthogroups, 'r'):
        for gene in genesids:
            if len(re.findall(gene+'.', l)) > 0:
                l2 = l.strip().split()
                for g in l2[1:]:
                    if g[:2] != 'AT':
                        ortho[gene].append(g.replace('mRNA', 'Gene').rsplit('.', maxsplit=1)[0])
    for k, v in ortho.items():
        ortho[k] = ','.join(set(v))
    ortho = pd.DataFrame({'athal_gene': ortho.keys(), 'cur_genes': ortho.values()})
    outgene = genes.merge(ortho, how='outer', on=['athal_gene'])
    outgene.to_csv(f'{cwd}/kk_orthology_marker_genes.txt', header=True, sep='\t', index=False)
    # </editor-fold>


    # <editor-fold desc="Get orthologs for Leaf atlas marker genes">
    # Markers from https://doi.org/10.1093/plcell/koaa060
    genesdf = pd.read_table(lagenefin, skiprows=2)
    genesdf.drop_duplicates(keep=False, inplace=True)
    genes = deque()
    cell = ''
    for row in genesdf.itertuples(index=False):
        if row[0][:2] == 'AT':
            genes.append([cell, row[0]])
        else:
            cell = row[0]
    genes = pd.DataFrame(genes)
    genes.columns = 'cell_type athal_gene'.split()
    genesids = set(genes.athal_gene)
    ortho = defaultdict(deque)
    for l in open(orthogroups, 'r'):
        for gene in genesids:
            if len(re.findall(gene+'.', l)) > 0:
                l2 = l.strip().split()
                for g in l2[1:]:
                    if g[:2] != 'AT':
                        ortho[gene].append(g.replace('mRNA', 'Gene').rsplit('.', maxsplit=1)[0])
    for k, v in ortho.items():
        ortho[k] = list(set(v))
    ortho = pd.DataFrame({'athal_gene': ortho.keys(), 'cur_genes': ortho.values()})
    outgene = genes.merge(ortho, how='outer', on=['athal_gene'])
    outgene = outgene.explode('cur_genes')
    outgene.to_csv(f'{cwd}/leafatlas_orthology_marker_genes.txt', header=True, sep='\t', index=False)
    # </editor-fold>


    # <editor-fold desc="Get orthologs for Leaf atlas marker genes 2">
    # Markers from https://doi.org/10.1093/plphys/kiab489
    genesdf = pd.read_csv(la2genefin)
    population = genesdf.Population
    for i in range(len(population)):
        if population[i] is np.nan:
            population[i] = population[i-1]
    genesdf.Population = population
    genes = genesdf['Population TAIR'.split()].copy()
    genes.columns = 'cell_type athal_gene'.split()
    genesids = set(genes.athal_gene)
    ortho = defaultdict(deque)
    for l in open(orthogroups, 'r'):
        for gene in genesids:
            if len(re.findall(gene+'.', l)) > 0:
                l2 = l.strip().split()
                for g in l2[1:]:
                    if g[:2] != 'AT':
                        ortho[gene].append(g.replace('mRNA', 'Gene').rsplit('.', maxsplit=1)[0])
    for k, v in ortho.items():
        ortho[k] = list(set(v))
    ortho = pd.DataFrame({'athal_gene': ortho.keys(), 'cur_genes': ortho.values()})
    outgene = genes.merge(ortho, how='outer', on=['athal_gene'])
    outgene = outgene.explode('cur_genes')
    outgene.to_csv(f'{cwd}/leafatlas2_orthology_marker_genes.txt', header=True, sep='\t', index=False)
    # </editor-fold>


    # <editor-fold desc="Get orthologs for three layers defined by Leaf atlas marker genes 2">
    # Markers from https://doi.org/10.1093/plphys/kiab489
    genes = pd.read_excel(la2layerfin, skiprows=1)
    genes = genes.melt()
    genes = genes.loc[~pd.isna(genes.value)]
    genes.value = [v.split('-')[0] for v in genes.value]
    genes.columns = 'cell_type athal_gene'.split()
    genesids = set(genes.athal_gene)
    ortho = defaultdict(deque)
    for l in open(orthogroups, 'r'):
        for gene in genesids:
            if len(re.findall(gene+'.', l)) > 0:
                l2 = l.strip().split()
                for g in l2[1:]:
                    if g[:2] != 'AT':
                        ortho[gene].append(g.replace('mRNA', 'Gene').rsplit('.', maxsplit=1)[0])
    for k, v in ortho.items():
        ortho[k] = list(set(v))
    ortho = pd.DataFrame({'athal_gene': ortho.keys(), 'cur_genes': ortho.values()})
    outgene = genes.merge(ortho, how='outer', on=['athal_gene'])
    outgene = outgene.explode('cur_genes')
    outgene.to_csv(f'{cwd}/leafatlas2_layer_orthology_marker_genes.txt', header=True, sep='\t', index=False)
    # </editor-fold>


    # <editor-fold desc="Get orthologs for adaxial and abaxial epidermis defined by Leaf atlas marker genes 2">
    # Markers from https://doi.org/10.1093/plphys/kiab489
    genes = pd.read_excel(la2epidermfin, skiprows=1)
    genes = genes.melt()
    genes = genes.loc[~pd.isna(genes.value)]
    genes.value = [v.split('-')[0] for v in genes.value]
    genes.columns = 'cell_type athal_gene'.split()
    genesids = set(genes.athal_gene)
    ortho = defaultdict(deque)
    for l in open(orthogroups, 'r'):
        for gene in genesids:
            if len(re.findall(gene+'.', l)) > 0:
                l2 = l.strip().split()
                for g in l2[1:]:
                    if g[:2] != 'AT':
                        ortho[gene].append(g.replace('mRNA', 'Gene').rsplit('.', maxsplit=1)[0])
    for k, v in ortho.items():
        ortho[k] = list(set(v))
    ortho = pd.DataFrame({'athal_gene': ortho.keys(), 'cur_genes': ortho.values()})
    outgene = genes.merge(ortho, how='outer', on=['athal_gene'])
    outgene = outgene.explode('cur_genes')
    outgene.to_csv(f'{cwd}/leafatlas2_epiderm_orthology_marker_genes.txt', header=True, sep='\t', index=False)
    # </editor-fold>


    # <editor-fold desc="Get structurally orthologous genes between Athal and Currot">
    ## Reads the output of foldseek structural alignments between Currot and Athal proteins
    ## and then finds the best structural orthologs athal gene for Currot genes
    cwd2 = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/structurome/'

    # Read uniprot to TAIR mapping
    uni2tair = pd.read_table(uniprot2tairfin, header=None)
    uni2tair[0] = uni2tair[0].str.replace("-.", '')
    uni2tair.drop_duplicates(inplace=True)

    cur2athal = defaultdict(deque)
    athal2cur = defaultdict(deque)
    with open(curmrnafin, 'r') as fin:
        for line in fin:
            cur2athal[line.strip()]
    for g in set(uni2tair[2]):
        athal2cur[g]

    alluni = pd.read_table(allathaluni)
    alluni = alluni[['Entry', 'Gene Names']]
    alluni.dropna(inplace=True)
    athalgenes = defaultdict(deque)
    for row in alluni.itertuples(index=False):
        for gene in row[1].replace('/', ' ').split():
            if gene[:2].upper() == 'AT' and len(gene) == 9:
                athalgenes[row[0]].append(gene.upper())
    for k in athalgenes:
        athalgenes[k] = set(athalgenes[k])

    # Add reciprocal best hits (RBH) between currot and athal
    df = pd.read_table(f'{cwd2}/athal_cur_rbh.tsv')
    df['query'] = df['query'].str.replace('.pdb.gz', '')
    df['target'] = df['target'].str.replace('.pdb.gz', '')
    df['target'] = df['target'].str.replace('-F1-model_v4', '')
    df['target'] = df['target'].str.replace('AF-', '')
    obsleleteuniids = {'B3H5K3': 'At5g47690'.upper(),
                       'Q9SN96': 'At3g46430'.upper(),
                       'Q9FFQ8': 'At5g63550'.upper(),
                       'A0A1P8ARR1': 'At1g80810'.upper(),
                       'F4KJ46': 'At5g53380'.upper(),
                       'Q9LZG3': 'At3g43850'.upper(),
                       'F4JHC9': 'At4g00520'.upper()}
    for row in df.itertuples(index=False):
        athalids = athalgenes[row[1]] #.loc[uni2tair[0] == row[1]][2].to_list()
        if len(athalids) == 0:
            try:
                cur2athal[row[0]].append(obsleleteuniids[row[1]])
                athal2cur[obsleleteuniids[row[1]]].append(row[0])
            except KeyError:
                print(row)
        for athal in athalids:
            cur2athal[row[0]].append(athal)
            athal2cur[athal].append(row[0])

    # Save RBH ids as they will not be modified
    rbhbkp = deepcopy(cur2athal)

    rbhcur = set([k for k in cur2athal if len(cur2athal[k]) > 0])
    rbhathal = set([k for k in athal2cur if len(athal2cur[k]) > 0])
    print(sum([len(v) for v in cur2athal.values()]), sum([len(v) for v in athal2cur.values()]))

    # Select best match from the alignment of cur mrna to the Athal
    df = pd.read_table(f'{cwd2}/athal_cur_search.tsv')
    df['query'] = df['query'].str.replace('.pdb.gz', '')
    df['target'] = df['target'].str.replace('.pdb.gz', '')
    df['target'] = df['target'].str.replace('-F1-model_v4', '')
    df['target'] = df['target'].str.replace('AF-', '')
    df = df.loc[~df['query'].isin(rbhcur)]
    for grp in df.groupby('query'):
        mid = grp[0]
        uniid = grp[1].iat[0,1]
        athalids = athalgenes[uniid]
        if len(athalids) == 0:
            try:
                cur2athal[mid].append(obsleleteuniids[uniid])
                if obsleleteuniids[uniid] not in rbhathal:
                    athal2cur[obsleleteuniids[uniid]].append(mid)
            except KeyError:
                print(mid)
        for athal in athalids:
            cur2athal[mid].append(athal)
            if athal not in rbhathal:
                athal2cur[athal].append(mid)

    # Select best match from the alignment of athal mrna to Cur
    df = pd.read_table(f'{cwd2}/cur_athal_search.tsv')
    df['target'] = df['target'].str.replace('.pdb.gz', '')
    df['query'] = df['query'].str.replace('.pdb.gz', '')
    df['query'] = df['query'].str.replace('-F1-model_v4', '')
    df['query'] = df['query'].str.replace('AF-', '')
    df = df.loc[~df['query'].isin(rbhathal)]
    for grp in df.groupby('query'):
        mid = grp[1].iat[0, 1]
        uniid = grp[0]
        athalids = athalgenes[uniid]
        if len(athalids) == 0:
            try:
                if mid not in rbhcur:
                    cur2athal[mid].append(obsleleteuniids[uniid])
                athal2cur[obsleleteuniids[uniid]].append(mid)
            except KeyError:
                print(uniid)
        for athal in athalids:
            if mid not in rbhcur:
                cur2athal[mid].append(athal)
            athal2cur[athal].append(mid)

    dfdata = deque()
    for k, vs in cur2athal.items():
        for v in vs:
            dfdata.append((k, v))
    for k, vs in athal2cur.items():
        for v in vs:
            dfdata.append((v, k))
    df = pd.DataFrame(dfdata).drop_duplicates()
    isrbh = deque()
    for row in df.itertuples(index=False):
        if row[1] in rbhbkp[row[0]]:
            isrbh.append(1)
        else:
            isrbh.append(0)
    df['rbh'] = isrbh
    df.sort_values(by=[0, 1], inplace=True)
    df.columns = 'currot_id athal_id reciprocal_match'.split()
    df.to_csv(f'{cwd2}/structural_orthologs_cur_athal.tsv', sep='\t', index=False)

    # df = pd.read_table(f'{cwd2}/structural_orthologs_cur_athal.tsv')
    df['currot_id'] = df['currot_id'].str.replace('mRNA', 'Gene')
    df['currot_id'] = df['currot_id'].str.replace(r'\..$', '')
    df.sort_values(by='reciprocal_match', inplace=True, ascending=False)
    df = df.drop_duplicates(subset='currot_id athal_id'.split())
    df.sort_values(by='currot_id athal_id'.split(), inplace=True)
    df.to_csv(f'{cwd2}/structural_orthologs_cur_athal.geneids.tsv', sep='\t', index=False)

    # </editor-fold>


    # <editor-fold desc="Get structural orthologs for KK marker genes">
    strortho = pd.read_table(strorthofin)
    genes = pd.read_table(kkgenefin, header=None, sep=',')
    genes.columns = 'cell_type athal_gene'.split()
    genesids = set(genes.athal_gene)
    df = strortho.loc[strortho.athal_id.isin(genesids)].copy()
    df = df.merge(genes, left_on="athal_id", right_on="athal_gene", how='outer')
    df = df['cell_type athal_gene currot_id reciprocal_match'.split()]
    df.sort_values(by='cell_type athal_gene currot_id reciprocal_match'.split(), inplace=True)
    df.to_csv(f'{cwd}/kk_structure_orthology_marker_genes.txt', index=False, sep='\t')
    # </editor-fold>


    # <editor-fold desc="Get structural orthologs for leaf-atlas marker genes">
    strortho = pd.read_table(strorthofin)
    genesdf = pd.read_table(lagenefin, skiprows=2)
    genesdf.drop_duplicates(keep=False, inplace=True)
    genes = deque()
    cell = ''
    for row in genesdf.itertuples(index=False):
        if row[0][:2] == 'AT':
            genes.append([cell, row[0]])
        else:
            cell = row[0]
    genes = pd.DataFrame(genes)
    genes.columns = 'cell_type athal_gene'.split()
    genesids = set(genes.athal_gene)
    df = strortho.loc[strortho.athal_id.isin(genesids)].copy()
    df = df.merge(genes, left_on="athal_id", right_on="athal_gene", how='outer')
    df = df['cell_type athal_gene currot_id reciprocal_match'.split()]
    df.sort_values(by='cell_type athal_gene currot_id reciprocal_match'.split(), inplace=True)
    df.to_csv(f'{cwd}/leafatlas_structure_orthology_marker_genes.txt', index=False, sep='\t')
    # </editor-fold>


    # <editor-fold desc="Get structural orthologs for leaf-atlas marker genes 2 ">
    strortho = pd.read_table(strorthofin)
    genesdf = pd.read_csv(la2genefin)
    population = genesdf.Population
    for i in range(len(population)):
        if population[i] is np.nan:
            population[i] = population[i-1]
    genesdf.Population = population
    genes = genesdf['Population TAIR'.split()].copy()
    genes.columns = 'cell_type athal_gene'.split()
    genesids = set(genes.athal_gene)
    df = strortho.loc[strortho.athal_id.isin(genesids)].copy()
    df = df.merge(genes, left_on="athal_id", right_on="athal_gene", how='outer')
    df = df['cell_type athal_gene currot_id reciprocal_match'.split()]
    df.sort_values(by='cell_type athal_gene currot_id reciprocal_match'.split(), inplace=True)
    df.to_csv(f'{cwd}/leafatlas2_structure_orthology_marker_genes.txt', index=False, sep='\t')
    # </editor-fold>\


    # <editor-fold desc="Get structural orthologs for leaf-atlas marker genes 2 ">
    strortho = pd.read_table(strorthofin)
    genesdf = pd.read_csv(la2genefin)
    population = genesdf.Population
    for i in range(len(population)):
        if population[i] is np.nan:
            population[i] = population[i-1]
    genesdf.Population = population
    genes = genesdf['Population TAIR'.split()].copy()
    genes.columns = 'cell_type athal_gene'.split()
    genesids = set(genes.athal_gene)
    df = strortho.loc[strortho.athal_id.isin(genesids)].copy()
    df = df.merge(genes, left_on="athal_id", right_on="athal_gene", how='outer')
    df = df['cell_type athal_gene currot_id reciprocal_match'.split()]
    df.sort_values(by='cell_type athal_gene currot_id reciprocal_match'.split(), inplace=True)
    df.to_csv(f'{cwd}/leafatlas2_structure_orthology_marker_genes.txt', index=False, sep='\t')
    # </editor-fold>\


    # <editor-fold desc="Get structural orthologs for three layers defined by Leaf atlas marker genes 2">
    strortho = pd.read_table(strorthofin)
    genes = pd.read_excel(la2layerfin, skiprows=1)
    genes = genes.melt()
    genes = genes.loc[~pd.isna(genes.value)]
    genes.value = [v.split('-')[0] for v in genes.value]
    genes.columns = 'cell_type athal_gene'.split()
    genesids = set(genes.athal_gene)
    df = strortho.loc[strortho.athal_id.isin(genesids)].copy()
    df = df.merge(genes, left_on="athal_id", right_on="athal_gene", how='outer')
    df = df['cell_type athal_gene currot_id reciprocal_match'.split()]
    df.sort_values(by='cell_type athal_gene currot_id reciprocal_match'.split(), inplace=True)
    df.to_csv(f'{cwd}/leafatlas2_layer_structure_orthology_marker_genes.txt', index=False, sep='\t')
    # </editor-fold>\


    # <editor-fold desc="Get structural orthologs for adaxial and abaxial epidermis defined by Leaf atlas marker genes 2">
    strortho = pd.read_table(strorthofin)
    genes = pd.read_excel(la2epidermfin, skiprows=1)
    genes = genes.melt()
    genes = genes.loc[~pd.isna(genes.value)]
    genes.value = [v.split('-')[0] for v in genes.value]
    genes.columns = 'cell_type athal_gene'.split()
    genesids = set(genes.athal_gene)
    df = strortho.loc[strortho.athal_id.isin(genesids)].copy()
    df = df.merge(genes, left_on="athal_id", right_on="athal_gene", how='outer')
    df = df['cell_type athal_gene currot_id reciprocal_match'.split()]
    df.sort_values(by='cell_type athal_gene currot_id reciprocal_match'.split(), inplace=True)
    df.to_csv(f'{cwd}/leafatlas2_epiderm_structure_orthology_marker_genes.txt', index=False, sep='\t')
    # </editor-fold>


    # <editor-fold desc="Get sequence and structural orthologous for thaliana DNA repair genes">
    dna_repair = pd.read_csv(f'TAIR_DNA_repair_gene_ids.txt', header=None)
    genesids = set(dna_repair[0])
    ortho = defaultdict(deque)
    for l in open(orthogroups, 'r'):
        for gene in genesids:
            if len(re.findall(gene+'.', l)) > 0:
                l2 = l.strip().split()
                for g in l2[1:]:
                    if g[:2] != 'AT':
                        ortho[gene].append(g.replace('mRNA', 'Gene').rsplit('.', maxsplit=1)[0])

    strortho = pd.read_table(strorthofin)
    df = strortho.loc[strortho.athal_id.isin(genesids)].copy()
    comb = set([(a, c) for a in genesids for c in ortho[a]] + [(a, c) for a in genesids for c in list(df.loc[df.athal_id == a, 'currot_id'])])
    comb = pd.DataFrame(comb)
    comb.sort_values([1], inplace=True)
    comb.to_csv('dna_repair_genes_currot.tsv', sep='\t', header=False, index=False)

    # </editor-fold>

    return
# END


def select_pdb_from_af2():
    """
    For each mRNA, select the highest scoring PDB (ranked_0) and save into a .tar file
    """
    import tarfile
    from tqdm import tqdm
    import os
    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/protein_structure/currot_pdbs/'
    with tarfile.open(f'{cwd}/selected_pdbs.tar', "w") as fout:
        mrnas = os.listdir(f'{cwd}/pdbs/')
        for mrna in tqdm(mrnas):
            try:
                fout.add(f'{cwd}/pdbs/{mrna}/ranked_0.pdb', arcname=f'{mrna}.pdb') # ranked 0 corresponds to the PDB with the highest score.
            except FileNotFoundError:
                pass
    return
# END

