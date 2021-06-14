#!/usr/bin/env python

CWD = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/rtiger_out/'
SAMPLES = ('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15')
CHRS = ('CUR1G', 'CUR2G', 'CUR3G', 'CUR4G', 'CUR5G', 'CUR6G', 'CUR7G', 'CUR8G')
INDIR = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/'
SNPFIN = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/strict_syn_snp.selected.txt'
CHRLEN = {'CUR1G': 46975282,
          'CUR2G': 33806098,
          'CUR3G': 26861604,
          'CUR4G': 26096899,
          'CUR5G': 18585576,
          'CUR6G': 27136638,
          'CUR7G': 25539660,
          'CUR8G': 23045982}

import os
import re
import pandas as pd
from collections import deque, defaultdict
from matplotlib import pyplot as plt
import numpy as np

def mergepdf(fins, fout):
    """
    Merge given PDF files in fins and save the combined file in fout.
    """
    from PyPDF2 import PdfFileMerger, PdfFileReader
    # Call the PdfFileMerger
    mergedObject = PdfFileMerger()
    for fin in fins: mergedObject.append(PdfFileReader(fin, 'rb'))
    # Write all the files into a file which is named as shown below
    mergedObject.write(fout)

################################################################################
df = pd.DataFrame()
for chr in CHRS:
    print(chr)
    os.chdir(CWD + chr + '/rtiger_co_q40/')     # Analyse samples with q=40
    # os.chdir(CWD + chr + '/rtiger_co/')       # Analyse samples with q=10
    dirs = [d for d in os.listdir() if os.path.isdir(d)]
    for d in dirs:
        s = [s for s in SAMPLES if s+'_' in d][0]
        # print(s)
        bc = re.findall(r'[a-zA-Z]{16}', d)[0]
        g = [f for f in os.listdir(d) if 'CompleteBlock-state' in f][0]
        bcdf = pd.read_table(d + '/' + g, header=None)
        bcdf.columns = ['chr', 'start', 'end', 'genotype']
        bcdf['sample'] = s
        bcdf['bc'] = bc
        bcdf = bcdf[['sample', 'bc', 'chr', 'start', 'end', 'genotype']]
        df = pd.concat([df, bcdf])


# Get the coordinate of the terminal (101st SNP) marker. COs identified before this marker would be filtered out.
term = pd.DataFrame()
for chr in CHRS:
    print(chr)
    dirpath = CWD + chr + "/input/"
    fins = os.listdir(dirpath)
    chrdata = deque()
    for f in fins:
        s = [s for s in SAMPLES if s+'_' in f][0]
        bc = re.findall(r'[a-zA-Z]{16}', f)[0]
        cnt = 0
        with open(dirpath + f, 'r') as fin:
            for line in fin:
                cnt += 1
                if cnt < 101: continue
                else:
                    chrdata.append([s, bc, chr, int(line.strip().split()[1])])
                    break
    term = pd.concat([term, pd.DataFrame(chrdata)])
term.columns = ['sample', 'bc', 'chr', 'cutoff']

SNPS = pd.read_table(SNPFIN, header=None)
SNPS = SNPS.drop([2], axis=1)
SNPS.columns = ['chr', 'pos', 'ref', 'alt']

########## Plot CO Counts in BCs ##################

cocnt = deque()
for grp in df.groupby(['sample', 'bc', 'chr']):
    # break
    if grp[1].shape[0] == 1:
        cocnt.append(list(grp[0]) + [0])
    else:
        garb = grp[1].copy()
        garb['tmp'] = (garb['end'][:-1].reset_index(drop=True) + garb['start'][1:].reset_index(drop=True))/2
        tmpdf = garb.merge(term, on=['sample', 'bc', 'chr']).copy()
        cocnt.append(list(grp[0]) + [sum(tmpdf['tmp'] > tmpdf['cutoff'])])

cocnt = pd.DataFrame(cocnt)
cocnt.columns = ['sample', 'bc', 'chr', 'co_count']
fig = plt.figure(figsize=[8, 10])
plti = 1
M = max(cocnt.co_count)
for grp in cocnt.groupby(['sample']):
    ax = fig.add_subplot(4, 1, plti)
    ax.set_ylim([0, 525])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    x_start = [0] * len(CHRS)
    for i in range(0, (M+1)):
        tmp = grp[1].loc[grp[1].co_count == i].copy()
        cnt = deque()
        for chr in CHRS:
            # cnt.append(garb.loc[garb.chr == chr].shape[0])
            cnt.append(tmp.loc[tmp.chr == chr].shape[0])
        ax.bar(CHRS, cnt, label=str(i), bottom=x_start, alpha=0.75)
        x_start = [x_start[j] + cnt[j] for j in range(len(cnt))]
    ax.set_xlabel(grp[0])
    ax.set_ylabel("Number of BCs")
    # ax.legend()
    plti += 1
plt.legend(loc=[0.1, -0.85], ncol=7, title="Numbe of COs")
plt.tight_layout()
plt.savefig(CWD + "co_counts_in_chromosomes_q40.pdf")


######## Plot CO Frequency along chromosome ################

coloc = pd.DataFrame()
for grp in df.groupby(['sample', 'bc', 'chr']):
    if grp[1].shape[0] == 1: continue
    loc = (grp[1]['end'][:-1].reset_index(drop=True) + grp[1]['start'][1:].reset_index(drop=True))/2
    loc = pd.DataFrame(loc)
    loc.columns = ['pos']
    loc['sample'] = grp[0][0]
    loc['bc'] = grp[0][1]
    loc['chr'] = grp[0][2]
    loc = loc[['sample', 'bc', 'chr', 'pos']]
    coloc = pd.concat([coloc, loc])

coloc2 = coloc.merge(term, on=['sample', 'bc', 'chr'])
coloc2 = coloc2.loc[coloc2['pos'] >= coloc2['cutoff']]
fig = plt.figure(figsize=[12, 9])
plti = 1
linestyle = {'WT_1': '-k', 'WT_19': '-.k', 'MUT_11_1': '-r', 'MUT_15': '-.r'}
for grp in coloc2.groupby('chr'):
    ax = fig.add_subplot(4, 2, plti)
    ax.set_ylim([0, 25])
    for grp2 in grp[1].groupby('sample'):
        ax.set_xlim([0, CHRLEN['CUR1G']])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        bins = list(range(1, CHRLEN[grp[0]], 1000000)) + [CHRLEN[grp[0]]]
        nbin = len(bins) - 1
        poscnt = pd.cut(sorted(grp2[1]['pos']), bins).value_counts().values
        x = (np.array(bins)[0:-1] + np.array(bins)[1:])/2
        ax.plot(x, poscnt, linestyle[grp2[0]], label=grp2[0])
        # ax.bar(x, poscnt, label=grp2[0], width=500000)
        # ax.text(0.7, 0.7, grp2[0], transform=ax.transAxes)
        ax.set_ylabel('#CO')
    ax.set_xlabel(grp[0])
    ax.legend(fontsize='x-small', loc=1)
    plti += 1
plt.tight_layout()
plt.savefig(CWD +'co_counts_along_chromosome_q40.pdf')


######## Plot Chromosome CO-Count vs Average read depth #################
read_depth = pd.DataFrame()
for sample in SAMPLES:
    sdf = pd.DataFrame()
    bcs = os.listdir("{}/{}/barcodes/".format(INDIR, sample))
    for bc in bcs:
        bcdf = deque()
        with open("{}/{}/barcodes/{}/mean_read_cov_q30_Q40.txt".format(INDIR, sample, bc), 'r') as fin:
            for line in fin:
                if 'CUR' in line:
                    line = line.strip().split()
                    bcdf.append([sample, bc, line[0], round(float(line[1]), 2)])
        bcdf = pd.DataFrame(bcdf)
        sdf = pd.concat([sdf, bcdf])
    read_depth = pd.concat([read_depth, sdf])
read_depth.columns = ['sample', 'bc', 'chr', 'depth']

cord = cocnt.merge(read_depth, on=['sample', 'bc', 'chr'])
import seaborn as sns
sns.catplot(x='co_count', y='depth', data=cord, hue='sample',
            col='chr', col_wrap=2, kind='violin', inner=None,
            height=3, aspect=2.5, linewidth=0)
plt.tight_layout()
plt.tight_layout()
plt.savefig(CWD+'CO_read_depth_distribution_q40.pdf')

################################################################################
######################## Find co-occuring COs ##################################
################################################################################
from matplotlib.ticker import AutoMinorLocator

hiloc = deque()
for grp in coloc2.groupby('chr'):
    tmp = grp[1].copy()
    pos = tmp['pos'].to_list()
    for p in pos:
        s = sum(abs(tmp['pos'] - p) < 20000)
        if s >= 1:
            hiloc.append([grp[0], p, s])
hiloc = pd.DataFrame(hiloc)
hiloc.columns = ['chr', 'pos', 'co']
fig = plt.figure(figsize=[25, 20])
plti = 1
for grp in hiloc.groupby('chr'):
    ax = fig.add_subplot(8, 1, plti)
    ax.set_ylim([0, 25])
    ax.set_xlim([0, CHRLEN['CUR1G']])
    ax.set_title(grp[0])
    ax.set_axisbelow(True)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(which='both', axis='both')
    ax.vlines(x=grp[1]['pos'], ymin=0, ymax=grp[1]['co'], linewidth=0.5, color='grey',zorder=1)
    ax.scatter(grp[1]['pos'], grp[1]['co'], color='black', s=4, zorder=3)
    ax.set_ylabel("CO Count")

    snps = SNPS.loc[SNPS['chr'] == grp[0]]
    ax.vlines(x=snps['pos'], ymin=0, ymax=25, color='lightblue', linewidth=0.1,zorder=0, alpha=0.1)
    plti+=1
    # break
ax.set_xlabel("Chromosome position (in Mbp)")
plt.tight_layout()
plt.tight_layout()
plt.savefig(CWD+'CO_clusters_q40.pdf')
plt.savefig(CWD+'CO_clusters_q40.png')


## Select plots for CO Clusters

hiloc2 = hiloc.loc[hiloc['co']>5].copy()
hiloc2.sort_values(['chr', 'pos'], inplace=True)
hiloc2.reset_index(inplace=True, drop=True)
posloc=[41]

c = hiloc2.iat[41, 0]
p = hiloc2.iat[41, 1]
fs = coloc2.loc[(coloc2['pos'] >= p-20000) & (coloc2['pos'] < p+20000) & (coloc2['chr']==c)]
fins = ['{CWD}/{c}/rtiger_co/{c}_{s}_{bc}_b30_q10.bt2.txt/GenotypePlot_{c}_{s}_{bc}_b30_q10.bt2.txt.pdf'.format(CWD=CWD, c=c, s=row[0], bc=row[1]) for row in fs.itertuples(index=False)]
fout = '/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/tests/temp.pdf'
mergepdf(fins, fout)



#############################################################################
############### Chromosomes with large Homozygous regions ###################
#############################################################################


def gethighhom(df, CWD, hetrat=0.5):
    fins = deque()
    cnt = 0
    for grp in df.groupby(['sample', 'bc', 'chr']):
        cnt += 1
        tl = sum(grp[1]['end'] - grp[1]['start'] + 1)
        het = grp[1].loc[grp[1]['genotype']=='AB']
        hetl = sum(het['end'] - het['start'] + 1)
        if hetl/tl < hetrat:
            fins.append('{CWD}/{c}/rtiger_co_q40/{c}_{s}_{bc}_b30_q40.bt2.txt/GenotypePlot_{c}_{s}_{bc}_b30_q40.bt2.txt.pdf'.format(CWD=CWD, c=grp[0][2], s=grp[0][0], bc=grp[0][1]))
    return (fins, cnt)

fins, cnt = gethighhom(df, CWD, 0.5)
# Number of Chromosomes tested: 11378 of which 526 (4.6%) have > 50% homozygous regions
mergepdf(fins, CWD+'high_homo_chrs_q40_hom50.pdf')

fins, cnt = gethighhom(df, CWD, 0.1)
# Number of Chromosomes tested: 11378 of which 168 (1.4%) have > 90% homozygous regions
mergepdf(fins, CWD+'high_homo_chrs_q40_hom90.pdf')
alls = [grp[0][0] for grp in df.groupby(['sample','bc', 'chr'])]
selecteds = [s for f in fins for s in SAMPLES if s+"_" in f]

allchrs = [grp[0][2] for grp in df.groupby(['sample','bc', 'chr'])]
selectedchrs = [c for f in fins for c in CHRS if c in f]

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(2, 2, 1)
ax.bar(SAMPLES, [Counter(alls)[k] for k in SAMPLES], label='Total')
ax.bar(SAMPLES, [Counter(selecteds)[k] for k in SAMPLES], label='Homo')
ax.set_ylabel('Chr Count')
ax.legend()
ax = fig.add_subplot(2, 2, 2)
ax.bar(SAMPLES, [Counter(selecteds)[k]*100/Counter(alls)[k] for k in SAMPLES], label='Total')
ax.set_ylabel('Percent Homo')

ax = fig.add_subplot(2, 2, 3)
ax.bar(CHRS, [Counter(allchrs)[k] for k in CHRS], label='Total')
ax.bar(CHRS, [Counter(selectedchrs)[k] for k in CHRS], label='Homo')
ax.set_ylabel('Chr Count')
ax.legend()
ax = fig.add_subplot(2, 2, 4)
ax.bar(CHRS, [Counter(selectedchrs)[k]*100/Counter(allchrs)[k] for k in CHRS], label='Total')
ax.set_ylabel('Percent Homo')
plt.tight_layout()
plt.savefig(CWD+'homo_dist_sample_chr.pdf')


################################################################################
############################## Plot CO Density #################################
################################################################################

fig = plt.figure(figsize=[8, 16])
plti = 1
for grp in coloc.groupby('chr'):
    ax = fig.add_subplot(9, 1, plti)
    data = sorted(grp[1]['pos'])
    bins = list(range(1, CHRLEN[grp[0]], 1000000)) + [CHRLEN[grp[0]]]
    y = pd.cut(data, bins).value_counts().values
    x = [(bins[i] + bins[i+1])/2 for i in range(len(bins) - 1)]
    ax.plot(x, y)

    # snps = SNPS.loc[SNPS['chr'] == grp[0]]
    # ax.vlines(x=snps['pos'], ymin=0, ymax=max(y), color='lightblue', linewidth=0.1,zorder=0, alpha=0.1)
    ax.set_xlim([0, CHRLEN['CUR1G']])
    ax.set_title(grp[0], size=10)
    ax.set_ylabel('#CO')
    plti+=1


# Read CUR6G reversed COs
os.chdir(CWD+'/CUR6G/reversed_rtiger_co_q40/')
dirs = [d for d in os.listdir() if os.path.isdir(d)]
revdf = pd.DataFrame()
for d in dirs:
    s = [s for s in SAMPLES if s+'_' in d][0]
    # print(s)
    bc = re.findall(r'[a-zA-Z]{16}', d)[0]
    g = [f for f in os.listdir(d) if 'CompleteBlock-state' in f][0]
    bcdf = pd.read_table(d + '/' + g, header=None)
    bcdf.columns = ['chr', 'start', 'end', 'genotype']
    bcdf['sample'] = s
    bcdf['bc'] = bc
    bcdf = bcdf[['sample', 'bc', 'chr', 'start', 'end', 'genotype']]
    revdf = pd.concat([revdf, bcdf])

revcoloc = pd.DataFrame()
for grp in revdf.groupby(['sample', 'bc', 'chr']):
    if grp[1].shape[0] == 1: continue
    loc = (grp[1]['end'][:-1].reset_index(drop=True) + grp[1]['start'][1:].reset_index(drop=True))/2
    loc = pd.DataFrame(loc)
    loc.columns = ['pos']
    loc['sample'] = grp[0][0]
    loc['bc'] = grp[0][1]
    loc['chr'] = grp[0][2]
    loc = loc[['sample', 'bc', 'chr', 'pos']]
    revcoloc = pd.concat([revcoloc, loc])

ax = fig.add_subplot(9, 1, plti)
data = sorted(revcoloc['pos'])
bins = list(range(1, CHRLEN['CUR6G'], 1000000)) + [CHRLEN['CUR6G']]
y = pd.cut(data, bins).value_counts().values
x = [(bins[i] + bins[i+1])/2 for i in range(len(bins) - 1)]
ax.plot(x, y)
ax.set_xlim([0, CHRLEN['CUR1G']])
ax.set_title('Reversed CUR6G', size=10)
ax.set_ylabel('#CO')
ax.set_xlabel('Chromosome Postion')
plt.tight_layout()

os.chdir(CWD)
plt.savefig('CO_density.pdf')


############################################################################
###### Check CO size to see if smaller COs can be filtered out #############
############################################################################

def mergegenotype(indf):
    outdf = pd.DataFrame()
    for grp in indf.groupby(['sample', 'bc', 'chr']):
        # if len(pd.unique(grp[1]['genotype'])) == 1: continue
        tmp = deque()
        current = ''
        for row in grp[1].itertuples(index=False):
        # for row in garb.itertuples(index=False):
            if current == '':
                current = list(row)
                continue
            if current[5] != row.genotype:
                tmp.append(current)
                current = list(row)
            else:
                current[4] = row.end
        tmp.append(current)
        tmp = pd.DataFrame(tmp)
        tmp.columns = grp[1].columns
        outdf = pd.concat([outdf, tmp])
    return outdf

df2 = pd.DataFrame()
for grp in df.groupby(['sample', 'bc', 'chr']):
    if grp[1].shape[0] == 1:
        if grp[1]['genotype'].to_list()[0] == 'AB':
            # print(grp[1])
            # break
            continue
    df2 = pd.concat([df2, grp[1]])

# Set all small CO regions as Het
df2['l'] = df2['end']- df2['start'] + 1
df2.loc[df2['l']<500000, 'genotype'] = 'AB'
df2 = mergegenotype(df2)
df2['l'] = df2['end']- df2['start'] + 1

# If small Het region is surrounded by same Hom region, then set Het -> Hom
df3 = pd.DataFrame()
for grp in df2.groupby(['sample', 'bc', 'chr']):
    tmp = grp[1].copy()
    if tmp.shape[0] == 1:
        df3 = pd.concat([df3, tmp])
        continue
    g = tmp['genotype'].to_numpy()
    heti = np.where(g=='AB')[0]
    s = tmp['l'].to_numpy()
    n = len(g)-1
    try:
        for i in heti:
            if s[i] > 500000: continue
            if i == 0: g[i] = g[i+1]
            elif i == n: g[i] = g[i-1]
            elif g[i-1] == g[i+1]:
                if g[i-1] == 'AB':
                    print('ERROR')
                g[i] = g[i-1]
    except:
        print(grp[1])
        break
    tmp['genotype'] = g
    df3 = pd.concat([df3, tmp])
df3 = mergegenotype(df3)
df3['l'] = df3['end']- df3['start'] + 1

garb = mergegenotype(df3)
garb['l'] = garb['end']- garb['start'] + 1




for row in df3.itertuples(index=False):
    break
sdf = pd.read_table('{CWD}{chr}/input_q40/{chr}_{sample}_{bc}_b30_q40.bt2.txt'.format(CWD=CWD, chr=row.chr, sample=row.sample, bc=row.bc), header=None)
sdf.columns = ['chr', 'pos', 'ref', 'ra', 'qry', 'qa']
garb = sdf.loc[(sdf['pos']>row.start) & (sdf['pos']<row.end)]














######## Plot CO Size Distribution #################

fig = plt.figure(figsize=[10, 8])
plti = 1
for grp in df.groupby(['sample']):
    ax = fig.add_subplot(4, 1, plti)
    for chr in CHRS:
        bad_chr = cocnt.loc[(cocnt['sample']==grp[0]) & (cocnt['chr']==chr) & (cocnt['co_count']==0)]
        tmp = grp[1].loc[grp[1]['chr'] == chr]
        tmp = tmp.loc[~(tmp['bc'].isin(bad_chr['bc']))]
        l = tmp['end'] - tmp['start'] + 1
        y, x, _ = ax.hist(l, bins=range(1, CHRLEN[chr], 500000), histtype='stepfilled', alpha=0, label=chr)
        ax.plot([j + 250000 for j in x[:-1]], y)
    ax.set_title(grp[0])
    ax.set_xlabel('Genotype Size')
    ax.set_ylabel("Frequency")

    # break
    plti += 1
plt.legend(loc=[0.1, -0.85], ncol=4, title="Chromosome")
plt.tight_layout()



################################################################################
############## Read coverage along chromosome for all samples ##################
################################################################################
snpsd = defaultdict(dict)
for row in SNPS.itertuples(index=False):
    snpsd[row[0]][row[1]] = 0

d = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/'
for sample in SAMPLES:
    ssnp = snpsd.copy()
    dirs = sorted([dir for dir in os.listdir(d + sample) if os.path.isdir(d+sample+'/'+dir)][:100])
    for bc in dirs:
        print(bc)
        with open(d+sample+'/'+bc+'/{sample}_{bc}_b30_q40.bt2.txt'.format(sample=sample, bc=bc), 'r') as f:
            for line in f:
                line = line.strip().split()
                ssnp[line[0]][int(line[1])] += int(line[3]) + int(line[5])

fig = plt.figure()
for chr in CHRS:
    ax = fig.add_subplot()
    pos = np.array(list(ssnp[chr].keys()), dtype='int')
    bins = list(range(1, CHRLEN[chr], 100000)) + [CHRLEN[chr]]
    bcnt = {}
    for i in range(len(bins) -1):
        b = bins[i]
        xs = np.argwhere((pos>b) & (pos < (b+100000)))
        s = 0
        for x in xs:
            s += ssnp[chr][int(pos[x])]
        try:
            bcnt[b] = s/(len(xs))
        except ZeroDivisionError as e:
            bcnt[b] = 0





    xbin = pd.cut(list(ssnp[chr].keys()), bins)
    x = sorted(np.random.choice(list(ssnp[chr].keys()), 10000))
    keys = sorted(list(bcnt.keys()))
    ax.plot(keys, [bcnt[k] for k in keys])































    snpsd = dict(zip(SNPS['chr'], SNPS['pos']))