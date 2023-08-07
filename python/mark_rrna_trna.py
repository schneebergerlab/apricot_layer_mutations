"""
Annotate rRNA and tRNA genes in the GFF files
"""
# <editor-fold desc="Define Import">
from collections import deque, defaultdict
import pandas as pd
from pybedtools import BedTool
from io import StringIO
# </editor-fold>

# <editor-fold desc="Define constants">
curgff = "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/EVM_PASA/pasa_on_mancur/cur.pasa_out.sort.gff3"
oragff = "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/ora/EVM_PASA/pasa_on_mancur/ora.pasa_out.sort.gff3"
rmfin = "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/repeat/RepeatMasker/cur.genome.v1.fasta.out"

# </editor-fold>

# <editor-fold desc="Supporting functions">
def readgff(f):
    # Read GFF data
    geneids = defaultdict(dict)
    gffdata = defaultdict(dict)
    first = True
    with open(f, 'r') as fin:
        for line in fin:
            if line.strip() == '': continue
            line = line.strip().split()
            if line[2] == 'gene':
                if first:
                    chr = line[0]
                    id = line[8].split(';')[0].split('=')[1]
                    geneids[line[0]][id] = int(line[3])
                    seq = '\t'.join(line) + '\n'
                    first = False
                else:
                    gffdata[chr][id] = pd.read_table(StringIO(seq), sep='\t', header=None)
                    chr = line[0]
                    id = line[8].split(';')[0].split('=')[1]
                    geneids[line[0]][id] = int(line[3])
                    seq = '\t'.join(line) + '\n'
            else:
                seq = seq + '\t'.join(line) + '\n'
        gffdata[chr][id] = pd.read_table(StringIO(seq), sep='\t', header=None)
    return geneids, gffdata
# END

# </editor-fold>

# <editor-fold desc="The following could be useful when for fixing annotations automatically">
# Read tRNA and rRNA position
annos = {'tRNA': deque(),
         'rRNA': deque()}
rnas = {'tRNA', 'rRNA'}
with open(rmfin, 'r') as rm:
    for i in range(3):
        rm.readline()
    for line in rm:
        line = line.strip().split()
        if line[10] in rnas:
            annos[line[10]].append((line[4], int(line[5]), int(line[6])))
annos['tRNA'] = pd.DataFrame(annos['tRNA'])
annos['rRNA'] = pd.DataFrame(annos['rRNA'])

correctedgff = defaultdict(dict)
for c, genes in gffdata.items():
    for gid, gdata in genes.items():
        cds = gdata.loc[gdata[2] == 'CDS', [0, 3, 4]]
        printdf(gdata)
        printdf(cds)
        break
    break

# </editor-fold>

# <editor-fold desc="Manual correction of rRNA genes">
currna = (('CUR3G', 10995306, 11584414), ('CUR4G', 21352409, 21491812))
geneids, gffdata = readgff(curgff)
with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/EVM_PASA/pasa_on_mancur/cur.pasa_out.sort.rrna.gff3', 'w') as fout:
    for chrom in sorted(list(geneids.keys())):
        for gid in sorted(geneids[chrom].keys(), key=lambda x: geneids[chrom][x]):
            if chrom not in {'CUR3G', 'CUR4G'}:
                fout.write(gffdata[chrom][gid].to_csv(header=False, index=False, sep='\t'))
            else:
                c, s, e = gffdata[chrom][gid].iloc[0, [0, 3, 4]]
                if chrom == 'CUR3G':
                    if e < 10995306 or s > 11584414:
                        fout.write(gffdata[chrom][gid].to_csv(header=False, index=False, sep='\t'))
                        continue
                if chrom == 'CUR4G':
                    if e < 21352409 or s > 21491812:
                        fout.write(gffdata[chrom][gid].to_csv(header=False, index=False, sep='\t'))
                        continue
                g = gffdata[chrom][gid].copy()
                if any(['protein_coding_gene' in i for i in g[8]]):
                    uniinfo = deque()
                    for i in g[8]:
                        if 'protein_coding_gene' in i:
                            uniinfo.append(f'{i.rsplit(";", 1)[0]};Name=rRNA')
                        else:
                            uniinfo.append(i)
                    g[8] = uniinfo
                    exon = g.loc[g[2] == 'exon'].iloc[0]
                    exon[8] = f'Parent={gid}'
                    g2 = g.iloc[0].copy()
                    g2[3: 5] = exon[3: 5]
                    fout.write(pd.DataFrame([g2, exon]).to_csv(header=False, index=False, sep='\t'))
                else:
                    fout.write(gffdata[chrom][gid].to_csv(header=False, index=False, sep='\t'))
                    continue


orarna = ('ORA3G', 13134842, 13706455)
geneids, gffdata = readgff(oragff)
with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/ora/EVM_PASA/pasa_on_mancur/ora.pasa_out.sort.rrna.gff3', 'w') as fout:
    for chrom in sorted(list(geneids.keys())):
        for gid in sorted(geneids[chrom].keys(), key=lambda x: geneids[chrom][x]):
            if chrom not in {'ORA3G'}:
                fout.write(gffdata[chrom][gid].to_csv(header=False, index=False, sep='\t'))
            else:
                c, s, e = gffdata[chrom][gid].iloc[0, [0, 3, 4]]
                if chrom == 'ORA3G':
                    if e < 13134842 or s > 13706455:
                        fout.write(gffdata[chrom][gid].to_csv(header=False, index=False, sep='\t'))
                        continue
                g = gffdata[chrom][gid].copy()
                if any(['protein_coding_gene' in i for i in g[8]]):
                    uniinfo = deque()
                    for i in g[8]:
                        if 'protein_coding_gene' in i:
                            uniinfo.append(f'{i.rsplit(";", 1)[0]};Name=rRNA')
                        else:
                            uniinfo.append(i)
                    g[8] = uniinfo
                    exon = g.loc[g[2] == 'exon'].iloc[0]
                    exon[8] = f'Parent={gid}'
                    g2 = g.iloc[0].copy()
                    g2[3: 5] = exon[3: 5]
                    fout.write(pd.DataFrame([g2, exon]).to_csv(header=False, index=False, sep='\t'))
                else:
                    fout.write(gffdata[chrom][gid].to_csv(header=False, index=False, sep='\t'))
                    continue



# </editor-fold>

