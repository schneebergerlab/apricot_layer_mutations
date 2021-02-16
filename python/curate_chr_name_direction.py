"""
Script to rename, reorder, and reorient the chromosomes to match conventions
"""
# Cur:
#   reverse: lg1, lg3, lg5
#   name:    lg1..8 -> CUR 6, 4, 3, 8, 2, 7, 5, 1

# Ora:
#   reverse: lg1, lg4, lg5, lg6, lg7, lg8
#   name:    lg1..8 -> ORA 6, 4, 3, 8, 2, 7, 5, 1

import os
import sys
print(sys.executable)
sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/')
from myUsefulFunctions import readfasta, writefasta
from Bio.Seq import reverse_complement
from collections import OrderedDict

lg_to_chr = {1: 6,
             2: 4,
             3: 3,
             4: 8,
             5: 2,
             6: 7,
             7: 5,
             8: 1}
# Correct CUR
to_rev = ['lg1', 'lg3', 'lg5']
os.chdir('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/manual_curated_assembly/cur/')
infasta = readfasta('cur.manual_curated.fasta')
outfasta = {}
for k,v in infasta.items():
    lg = k.split('_')[0]
    if lg in to_rev: outseq = reverse_complement(v)
    else: outseq = v
    outchr = 'CUR' + str(int(lg_to_chr[int(lg[-1])])) + 'G'
    outfasta[outchr] = outseq
outdict = OrderedDict()
sorted_keys = sorted(outfasta.keys())
for k in sorted_keys:
    outdict[k] = outfasta[k]
writefasta(outdict, 'cur.manual_curated.renamed.fasta')


# Correct ORA
to_rev = ['lg1', 'lg4', 'lg5', 'lg6', 'lg7', 'lg8']
os.chdir('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/manual_curated_assembly/ora/')
infasta = readfasta('ora.manual_curated.fasta')
outfasta = {}
for k,v in infasta.items():
    lg = k.split('_')[0]
    if lg in to_rev: outseq = reverse_complement(v)
    else: outseq = v
    outchr = 'ORA' + str(int(lg_to_chr[int(lg[-1])])) + 'G'
    outfasta[outchr] = outseq
outdict = OrderedDict()
sorted_keys = sorted(outfasta.keys())
for k in sorted_keys:
    outdict[k] = outfasta[k]
writefasta(outdict, 'ora.manual_curated.renamed.fasta')

