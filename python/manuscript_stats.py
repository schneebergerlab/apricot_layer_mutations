import sys
sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/software/hometools/')
from myUsefulFunctions import mergeRanges, subranges
import pandas as pd
from collections import deque

# Get strictly_syntenic and SR regions
syriout = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out', header=None)
synout = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syn_range.txt'
srout = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/sr_range.txt'
with open(synout, 'w') as syno, open(srout, 'w') as sro:
    for g in syriout.groupby(0):
        if g[0] == '-': continue
        syn = g[1].loc[syriout[10]=='SYN']
        syn[[1, 2]] = syn[[1, 2]].astype(int)
        sr = g[1].loc[syriout[10].isin(['INV', 'TRANS', 'INVTR', 'DUP', 'INVDP'])]
        sr[[1, 2]] = sr[[1, 2]].astype(int)
        synr = mergeRanges(np.array(list(zip(syn[1], syn[2]))))
        srr = mergeRanges(np.array(list(zip(sr[1], sr[2]))))
        strctsyn = subranges(synr, srr)
        for i in strctsyn:
            syno.write("{}\t{}\t{}\n".format(g[0], i[0], i[1]))
        for i in srr:
            sro.write("{}\t{}\t{}\n".format(g[0], i[0], i[1]))
