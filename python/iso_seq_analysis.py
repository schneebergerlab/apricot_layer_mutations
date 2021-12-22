import os
import pandas as pd
from collections import deque, Counter
import sys
sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/software/hometools/')
from myUsefulFunctions import revcomp
from matplotlib import pyplot as plt
import warnings
import seaborn as sns
import numpy as np
import pysam
import pathlib

SAMPLES=('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15')
INDIR='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/'
CUTS=(10, 25, 50, 100)

# Get BC read-count stats
for s in SAMPLES:
    bcrc = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/barcodes/{}_bc_readcnts.txt'.format(s), header=None, delimiter=' ')
    df = deque()
    with open(INDIR+s+"/cnt.txt", 'r') as fin:
        for line in fin:
            line = line.strip().split()
            df.append(line)
    df = pd.DataFrame(df)
    df[1] = df[1].str.replace('XC:Z:', '')
    df[1] = list(map(revcomp, df[1]))
    df[0] = df[0].astype(int)
    isoillrc = pd.merge(bcrc, df, how='inner', on=[1])
    fig=plt.figure(figsize=(8, 8))
    i = 1
    ax = fig.add_subplot(3,2,i)
    ax.scatter(isoillrc['0_x'], isoillrc['0_y'], s=1)
    ax.set_xlabel("scRNA read count")
    ax.set_ylabel("scIso-Seq read count")
    ax.set_title("Read count in cluster BC\nTotal Iso-Seq reads: {}".format(sum(isoillrc['0_y'])), fontdict={'fontsize': 8})
    i+=1
    ax = fig.add_subplot(3,2,i)
    ax.hist(isoillrc['0_y'], bins=100)
    ax.set_xlabel("scIso-Seq read count")
    ax.set_ylabel("Number of BC")
    i+=1
    for c in CUTS:
        g = df.loc[(df[0]<1000) & (df[0]>=c)]
        ax = fig.add_subplot(3, 2, i)
        ax.set_title(i)
        ax.hist(g[0], bins=100)
        ax.set_xlim(0, 1000)
        ax.set_title("Number of reads: {}\nNumber of cells: {}\ncutoff: {}".format(sum(g[0]), len(g[1]), c), fontdict={'fontsize':8})
        i+=1
    plt.tight_layout(h_pad=3, w_pad=3)
    plt.savefig(INDIR+s+"/bc_read_dist.pdf")
    plt.savefig(INDIR+s+"/bc_read_dist.png", transparent=True)
    plt.close()

# Filter BC reads and save into Fasta
for s in SAMPLES:
    bcclt = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/barcodes/{}_bcs_clstr_id.txt'.format(s), header=None, delimiter=' ')
    bcs = set(bcclt[1])
    bcsrev = set([revcomp(bc) for bc in bcs])
    clstrs = set(bcclt[0])
    bcdict = dict(zip(bcclt[1], bcclt[0]))
    outfa = {c: open(INDIR+'/{}/clstrs_{}_reads.fa'.format(s, c), 'w') for c in clstrs}
    b = pysam.AlignmentFile(INDIR+'/{}/{}.XC_sorted.bam'.format(s, s), mode='rb', check_sq=False)
    for read in b:
        if read.get_tag('XC') in bcsrev:
            outfa[bcdict[revcomp(read.get_tag('XC'))]].write(">{};{}\n{}\n".format(read.query_name, revcomp(read.get_tag('XC')), read.seq))
# with open(INDIR+'/{}/{}_reads.fa'.format(s, s), 'w') as fout:
    for v in outfa.values():
        v.close()

# Plot Alt-readcount frequency
CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/'
snps=pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mutations.regions', header=None)
snps = snps.head(39)
snps[1] = snps[1].astype(int)
snps.sort_values([7, 0, 1], inplace=True)
snps['af'] = ['high' if af >0.2 else 'low' for af in snps[6]]
BASE_DICT = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
snprc = pd.DataFrame()
for s in ('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15'):
    clsrc = sorted([i for i in os.listdir('{}/{}'.format(CWD, s)) if 'rc.txt' in i])
    for cls in clsrc:
        df = pd.read_table('{}/{}/{}'.format(CWD, s, cls), header=None)
        df['s'] = s
        df['cls'] = cls.split('.')[0]
        snprc = pd.concat([snprc, df])
snprc = snprc.loc[(snprc[0].astype(str)+snprc[1].astype(str)).isin(list(snps[0].astype(str)+snps[1].astype(str)))]

snpclsrc = deque()
snpclsac = deque()
snpclsaf = deque()

warnings.filterwarnings("error")
for g in snprc.groupby(['s', 'cls']):
    b=0
    clsrc = deque()
    clsac = deque()
    clsaf = deque()
    for snp in snps.itertuples(index=False):
        snpg = g[1].loc[(g[1][0] == snp[0]) & (g[1][1] == snp[1])]
        if snpg.shape[0] == 0:
            clsrc.append(-1)
            clsac.append(-1)
            clsaf.append(0)
            continue
        if snpg.shape[0] != 1:
            print("ERROR in snp and cls: {} {} {}".format(snp[0]+str(snp[1]), g[0][0], g[0][1]))
            continue
        if snpg.iat[0, 3] == 0:
            clsrc.append(0)
            clsac.append(0)
            clsaf.append(0)
        else:
            clsrc.append(snpg.iat[0, 3])
            clsac.append(snpg.iat[0, BASE_DICT[snp[4]]])
            clsaf.append(snpg.iat[0, BASE_DICT[snp[4]]]/snpg.iat[0, 3])
    snpclsrc.append(g[0] + tuple(clsrc))
    snpclsac.append(g[0] + tuple(clsac))
    snpclsaf.append(g[0] + tuple(clsaf))

rcdf = pd.DataFrame(snpclsrc)
acdf = pd.DataFrame(snpclsac)
afdf = pd.DataFrame(snpclsaf)
xticks = list(snps.iloc[:, 7].astype(str) + "_" + snps.iloc[:, 0].astype(str)+":"+snps.iloc[:, 1].astype(str) + "_" + snps['af'])
yticks = list(rcdf.iloc[:, 0].astype(str) + "_" + rcdf.iloc[:, 1].astype(str))
fig = plt.figure(figsize=[10, 14])
plt.rcParams.update({'font.size': 8})
ax = fig.add_subplot(3, 1, 1)
ax = sns.heatmap(rcdf.iloc[:, 2:], linewidths=0.2, linecolor='k', center=0, cmap='seismic', xticklabels=False, yticklabels=yticks, cbar_kws={'fraction': 0.05})
ax.set_title("Read coverage")
ax = fig.add_subplot(3, 1, 2)
ax = sns.heatmap(acdf.iloc[:, 2:], linewidths=0.2, linecolor='k', center=0, cmap='seismic', xticklabels=False, yticklabels=yticks, cbar_kws={'fraction': 0.05})
ax.set_title("Alt allele read count")
ax = fig.add_subplot(3, 1, 3)
ax = sns.heatmap(afdf.iloc[:, 2:], linewidths=0.2, linecolor='k', center=0, cmap='seismic', xticklabels=xticks, yticklabels=yticks, cbar_kws={'fraction': 0.05})
ax.set_title("Alt allele frequency")
plt.tight_layout()
plt.savefig(CWD+"mut_rc_iso_seq.pdf")
plt.savefig(CWD+"mut_rc_iso_seq.png")
plt.close()

# Get cluster-specific read stats
s_cls_rc = {}
for s in SAMPLES:
    s_cls_rc[s] = {}
    # bams = [os.path.basename(i) for i in pathlib.Path(CWD+s).glob("clstrs*bam")]
    bams = [i for i in pathlib.Path(CWD+s).glob("clstrs*bam")]
    for bam in bams:
        with pysam.AlignmentFile(bam, "rb" ) as b:
            bname = os.path.basename(bam)
            rsize = deque()
            flags = deque()
            rname = deque()
            for r in b:
                rsize.append(r.infer_read_length())
                flags.append(r.flag)
                rn = r.reference_name if r.reference_name is not None else '*'
                rname.append(rn)
        s_cls_rc[s][bname.split('.')[0]] = {'rsize': tuple(rsize),
                                            'flags': tuple(flags),
                                            'rname': tuple(rname)}
        print(bname)

## Flags used: 0, 4, 16
for s, v in s_cls_rc.items():
    for c, v1 in v.items():
        print(s, c, Counter(v1[1]))



genecnts = Counter(pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.pasa_out.sort.protein_coding.3utr.bed', header=None)[0])
fig = plt.figure(figsize=[10,10])
i = 0
for s in SAMPLES:
    df = pd.DataFrame()
    for c, v in s_cls_rc[s].items():
        garb = pd.DataFrame(v)
        garb['clstr'] = c
        df = pd.concat([df, garb])
    df.sort_values(['clstr', 'rname'], inplace=True)
    df2 = df.loc[df.rname.str.contains('CUR')]
    rncnt = Counter(df2.rname)

    ax1 = fig.add_subplot(4,5,i+1)
    ax1 = sns.barplot(x=list(rncnt.keys()),y=list(rncnt.values()))
    ax1.set_ylabel("Read count")
    ax1.set_ylim([0,200000])
    ax1.tick_params(axis='x', labelrotation = 90)
    ax1.ticklabel_format(axis='y', style='scientific', scilimits=(0,2))

    ax2 = fig.add_subplot(4,5,i+2)
    ax2 = sns.barplot(x=list(rncnt.keys()),y=[rncnt[k]/genecnts[k] for k in rncnt.keys()])
    ax2.set_ylim([0,30])
    ax2.set_ylabel("Read per gene")
    ax2.tick_params(axis='x', labelrotation = 90)

    ax3 = fig.add_subplot(4,5,i+3)
    ax3 = sns.countplot(x='clstr', data=df2)
    ax3.set_ylim([0,200000])
    ax3.set_ylabel("Read count")
    ax3.set_xlabel(None)
    ax3.tick_params(axis='x', labelrotation = 90)
    ax3.ticklabel_format(axis='y', style='scientific', scilimits=(0,2))

    ax4 = fig.add_subplot(4,5,i+4)
    ax4.tick_params(axis='x', labelrotation = 90)
    ax4.set_ylim([0,6000])
    ax4 = sns.violinplot(x='clstr', y='rsize', data=df2, inner=None, color="white", cut=0)
    ax4.set_ylabel("Read size")
    ax4.set_xlabel(None)
    ax4.ticklabel_format(axis='y', style='scientific', scilimits=(0,2))

    ax5 = fig.add_subplot(4,5,i+5, frameon=False)
    ax5.axis('off')
    ax5.text(0,0,"{}\nTotal Reads:\n{}\nMapped Reads:\n{}".format(s, df.shape[0], sum(df.rname != '*')))
    i+=5

plt.tight_layout()
plt.savefig(CWD+"iso_read_mapping_stats.pdf")
plt.savefig(CWD+"iso_read_mapping_stats.png")
plt.close()