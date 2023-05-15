import os
import pandas as pd
from collections import deque, Counter, defaultdict
import sys
sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/software/hometools/')
# from myUsefulFunctions import revcomp
from hometools.hometools import revcomp
from matplotlib import pyplot as plt
import warnings
import seaborn as sns
import numpy as np
import pysam
import pathlib
from itertools import product
from matplotlib import colors as mcolors


def get_transcriptome_variants():
    """
    Compare the transcriptome for the four brancehs and check if there are genes
    with different transcripts/isoforms
    """
    import os
    import pandas as pd
    import pybedtools as bt
    from collections import defaultdict, deque
    from venn import venn
    import numpy as np
    import seaborn as sns

    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/allele_specific/variant_transcripts/'
    smfin = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.filtered.txt'
    gff = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/EVM_PASA/pasa_on_mancur/cur.pasa_out.sort.protein_coding.3utr.gff3'
    transdir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/allele_specific/'
    samples = 'WT_1 WT_19 MUT_11_1 MUT_15'.split()
    wsize = 10000

    # Change working directory
    os.chdir(cwd)

    # Get somatic mutations in bed format and then parse geneids near the somatic mutations
    sms = pd.read_table(smfin)
    df = pd.DataFrame({'chromosome': sms.chromosome, 'start': sms.position - 1, 'end': sms.position}).drop_duplicates()
    sms_bed = bt.BedTool.from_dataframe(df)
    muts = sms_bed.window(b=gff, w=wsize).saveas().filter(func=lambda x: x[5] == 'gene').saveas().to_dataframe()
    geneids = set(muts.apply(func=lambda x: x[11].split(';')[0].split('=')[1], axis=1).to_list())

    # Get transcripts that are matching geneids
    strid = {}  # Sample overlapping transcripts
    getdict = lambda x: {a.split()[0]: a.split()[1].replace('"', '') for a in x.split(';')[:-1]}
    for sample in samples:
        trids = sms_bed.window(b=f'{transdir}/{sample}/{sample}.resucue_rescued.gtf', w=wsize).saveas().filter(func=lambda x: x[5] == 'transcript').saveas().to_dataframe()
        # trids = set(trids.apply(func=lambda x: x[11].split(';')[1].strip().split()[1].replace('"', ''), axis=1).to_list())
        trids = set(trids.apply(func=lambda x: getdict(x[11])['transcript_id'], axis=1).to_list())
        tridsgene = {}
        with open(f'{transdir}/{sample}/{sample}.ML_MLresult_classification.txt', 'r') as fin:
            for line in fin:
                line = line.strip().split()
                if line[0] in trids:
                    tridsgene[line[0]] = [line[6], line[5], line[14]]
        # Transcripts and Genes for reference transcripts that are added in the
        # final GTF but absent in the classification file
        for t in trids:
            if t not in tridsgene:
                tridsgene[t] = t.replace('mRNA', 'Gene')
        strid[sample] = [trids, tridsgene]
    genetrs = {}
    for sample in samples:
        genetrs[sample] = defaultdict(deque)
        for k, v in strid[sample][1].items():
            if v[0] in geneids:
                genetrs[sample][v[0]].append(k)
    allgenes = set([k  for sample in samples for k in list(genetrs[sample].keys())])

    # Remove genes that only have full-splice match transcripts
    badgene = deque()
    for gene in allgenes:
        cat = deque()
        for sample in samples:
            try:
                for t in genetrs[sample][gene]:
                    cat.append(strid[sample][1][t][1])
            except KeyError:
                pass
        cat = set(cat)
        if cat == {'full-splice_match'}:
            badgene.append(gene)
    allgenes = allgenes.difference(set(badgene))

    # Remove genes that only have mono-exon
    badgene = deque()
    for gene in allgenes:
        cat = deque()
        for sample in samples:
            try:
                for t in genetrs[sample][gene]:
                    cat.append(strid[sample][1][t][2])
            except KeyError:
                pass
        cat = set(cat)
        if all([c in {'mono-exon', 'mono-exon_by_intron_retention', 'reference_match'} for c in cat]):
            badgene.append(gene)
    allgenes = allgenes.difference(set(badgene))

        for sample in samples:
            try:
                for t in genetrs[sample]['Gene.9100']:
                    print(sample, t, strid[sample][1][t])
            except KeyError:
                pass
    return
# END

def get_iso_seq_stats():
    """
    Compare the coverage/read-count (etc) of iso-seq reads to the scrna seq reads
    (Older code (before Anshupa's analysis and removal of duplicates by isoseq3)
     is overwritten).

    scrna read count and clusters generated using Anshupa's analysis
    """
    from hometools.hometools import revcomp
    sdict = {'WT_1': 'wt1',
             'WT_19': 'wt19',
             'MUT_11_1': 'mut_11_1',
             'MUT_15': 'mut_15'}
    indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/'
    CUTS = (10, 25, 50, 100)

    # Get BC read-count stats
    for k, v in sdict.items():
        bcrc = pd.read_table(f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/sahu_analysis/analysis/{v}_bc_readcnts.txt'.format(v), header=None, delimiter=' ')
        df = deque()
        with open(indir+k+"/cnt.txt", 'r') as fin:
            for line in fin:
                line = line.strip().split()
                df.append(line)
        df = pd.DataFrame(df)
        df[1] = df[1].str.replace('CB:Z:', '')
        df[1] = list(map(revcomp, df[1]))           # CB in Iso-Seq is reverse complemented
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



    return
# END
get_iso_seq_stats()

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
cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/'
# snps=pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mutations.regions', header=None)
muts = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.filtered.txt')
sdict = dict(zip(('WT_1', 'wt7', 'wt18', 'WT_19', 'MUT_11_1', 'mut11_2', 'mut4', 'MUT_15'), ('wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_4', 'mut_15')))
basedict = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
# snps = snps.head(39)
# snps[1] = snps[1].astype(int)
muts.sort_values(['branch', 'chromosome', 'position'], inplace=True)
pos = set(zip(muts.chromosome, muts.position, muts.alt_allele))
posdict = {}
for g in muts.groupby(['chromosome', 'position', 'alt_allele']):
    posdict[g[0]] = g[1].branch.to_list()
# snps['af'] = ['high' if af > 0.2 else 'low' for af in snps[6]]
mutsrc = pd.DataFrame()

for bname in ('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15'):
    clsrc = sorted([i for i in os.listdir('{}/{}'.format(cwd, bname)) if 'layer.rc.txt' in i])
    print(clsrc)
    bdf = pd.DataFrame()
    for cls in clsrc:
        clsdf = defaultdict(dict)
        clsname = cls.split(".")[0]
        with open('{}/{}/{}'.format(cwd, bname, cls), 'r') as fin:
            for line in fin:
                line = line.strip().split()
                if len(line) == 0: continue
                alts = deque()
                # print(line)
                for p in pos:
                    if (line[0], int(line[1])) == (p[0], p[1]):
                        alts.append(p[2])
                # print(alts)

                for alt in alts:
                    if int(line[3]) == 0:
                        clsdf[line[0], int(line[1]), alt] = {'rc': 0, 'ac': 0, 'af': 0}
                        continue
                    try:
                        clsdf[line[0], int(line[1]), alt] = {'rc': int(line[3]), 'ac': int(line[basedict[alt]]), 'af': round(int(line[basedict[alt]])/int(line[3]), 2)}
                    except KeyError:
                        try:
                            i = line.index(alt)
                            clsdf[line[0], int(line[1]), alt] = {'rc': int(line[3]), 'ac': int(line[i+1]), 'af': round(int(line[i+1])/int(line[3]), 2)}
                        except ValueError:
                            clsdf[line[0], int(line[1]), alt] = {'rc': int(line[3]), 'ac': 0, 'af': 0}
        clsdf = pd.DataFrame(clsdf).transpose()
        clsdf['cls'] = clsname
        clsdf['branch'] = sdict[bname]
        mutsrc = pd.concat([mutsrc, clsdf])
mutsrc.reset_index(level=[0, 1, 2], inplace=True)
mutsrc.columns = ['chromosome', 'position', 'alt_allele'] + list(mutsrc.columns[3:])
mutmerge = mutsrc.merge(muts, on=['chromosome', 'position', 'alt_allele', 'branch'], how='left')
mutmerge.sort_values(['chromosome', 'position', 'alt_allele'], inplace=True)
mutmerge.reset_index(drop=True, inplace=True)
mutmerge['p'] = mutmerge.chromosome.astype(str) + '_' + mutmerge.position.astype(str) + '_' + mutmerge.alt_allele.astype(str)
mutmerge['c'] = mutmerge.branch.astype(str) + '_' + mutmerge.cls.astype(str)
nonzeropos = set(mutmerge.loc[mutmerge.rc > 0].p)
mutmerge = mutmerge.loc[mutmerge.p.isin(nonzeropos)]


allpos = pd.DataFrame(product(mutmerge.p.unique(), mutmerge.c.unique()))
allpos.columns = ['p', 'c']
allpos = allpos.merge(mutmerge, how='left', on=['p', 'c'])
allpos.loc[allpos.rc.isna(), ['rc', 'ac', 'af']] = -1


L1pos = muts.loc[(muts.Layer=='L1') & (muts.branch.isin(['mut_11_1', 'mut_15', 'wt_1', 'wt_19']))]
L1pos = set(L1pos.chromosome.astype(str) + '_' + L1pos.position.astype(str) + '_' + L1pos.alt_allele.astype(str))
L2pos = muts.loc[(muts.Layer=='L2') & (muts.branch.isin(['mut_11_1', 'mut_15', 'wt_1', 'wt_19']))]
L2pos = set(L2pos.chromosome.astype(str) + '_' + L2pos.position.astype(str) + '_' + L2pos.alt_allele.astype(str))

cmap = [(1, 1, 1, 1) for i in np.linspace(0.1, 0, 10)] + [(1, 1-i, 1-i, 1) for i in np.linspace(0, 1, 110)]
cmap = mcolors.ListedColormap(cmap) # make color map

# TODO: consider adding sample information in the plot
fig = plt.figure(figsize=[8, 14])
ax = fig.add_subplot(3, 1, 1)
# Normalise values of rc, ac, af to be in range
rchm = allpos.pivot(index='c', columns='p')['rc']
rchm = pd.concat([rchm.loc[:, [c for c in rchm.columns if c in L1pos]], rchm.loc[:, [c for c in rchm.columns if c in L2pos]]], axis=1)
rchm = rchm*10/rchm.values.max()
rchm[rchm < 0] = -1
sns.heatmap(rchm, linewidths=0.1, linecolor='w', cmap=cmap, xticklabels=False, yticklabels=rchm.index, cbar_kws={'label': 'Normalized read count', 'fraction': 0.05}, ax=ax, vmin=-1)
ax.hlines([7, 14, 23], *ax.get_xlim(), color='k')
ax.vlines([len([c for c in rchm.columns if c in L1pos])], *ax.get_ylim(), color='k')
ax.set_ylabel('')
ax.set_xlabel('')
ax.set_title('Normalized read count')

ax = fig.add_subplot(3, 1, 2)
achm = allpos.pivot(index='c', columns='p')['ac']
achm = pd.concat([achm.loc[:, [c for c in achm.columns if c in L1pos]], achm.loc[:, [c for c in achm.columns if c in L2pos]]], axis=1)
achm = achm*10/achm.values.max()
achm[achm < 0] = -1
sns.heatmap(achm, linewidths=0.1, linecolor='w', cmap=cmap, xticklabels=False, yticklabels=achm.index, cbar_kws={'label': 'Normalized allele count', 'fraction': 0.05}, ax=ax, vmin=-1)
ax.hlines([7, 14, 23], *ax.get_xlim(), color='k')
ax.vlines([len([c for c in achm.columns if c in L1pos])], *ax.get_ylim(), color='k')
ax.set_ylabel('')
ax.set_xlabel('')
ax.set_title('Normalized allele count')

ax = fig.add_subplot(3, 1, 3)
afhm = allpos.pivot(index=['c'], columns='p')['af']
afhm = pd.concat([afhm.loc[:, [c for c in afhm.columns if c in L1pos]], afhm.loc[:, [c for c in afhm.columns if c in L2pos]]], axis=1)
afhm = afhm*10/afhm.values.max()
afhm[afhm < 0] = -1
sns.heatmap(afhm, linewidths=0.1, linecolor='w', cmap=cmap, xticklabels=True, yticklabels=afhm.index, cbar_kws={'label': 'Normalized allele frequency', 'fraction': 0.05}, ax=ax, vmin=-1)
# sns.heatmap(afhm, linewidths=0.1, linecolor='w', cmap=cmap, norm=norm, xticklabels=True, yticklabels=afhm.index, cbar_kws={'label': 'Normalized allele frequency', 'fraction': 0.05}, ax=ax, vmin=-1)
ax.hlines([7, 14, 23], *ax.get_xlim(), color='k')
ax.vlines([len([c for c in afhm.columns if c in L1pos])], *ax.get_ylim(), color='k')
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title('Normalized allele frequency')
plt.tight_layout()
plt.savefig('/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/var_af_in_iso_seq_clusters.pdf')
plt.close()

# snpclsrc = deque()
# snpclsac = deque()
# snpclsaf = deque()
#
# warnings.filterwarnings("error")
# for g in snprc.groupby(['s', 'cls']):
#     b=0
#     clsrc = deque()
#     clsac = deque()
#     clsaf = deque()
#     for snp in snps.itertuples(index=False):
#         snpg = g[1].loc[(g[1][0] == snp[0]) & (g[1][1] == snp[1])]
#         if snpg.shape[0] == 0:
#             clsrc.append(-1)
#             clsac.append(-1)
#             clsaf.append(0)
#             continue
#         if snpg.shape[0] != 1:
#             print("ERROR in snp and cls: {} {} {}".format(snp[0]+str(snp[1]), g[0][0], g[0][1]))
#             continue
#         if snpg.iat[0, 3] == 0:
#             clsrc.append(0)
#             clsac.append(0)
#             clsaf.append(0)
#         else:
#             clsrc.append(snpg.iat[0, 3])
#             clsac.append(snpg.iat[0, BASE_DICT[snp[4]]])
#             clsaf.append(snpg.iat[0, BASE_DICT[snp[4]]]/snpg.iat[0, 3])
#     snpclsrc.append(g[0] + tuple(clsrc))
#     snpclsac.append(g[0] + tuple(clsac))
#     snpclsaf.append(g[0] + tuple(clsaf))
#
# rcdf = pd.DataFrame(snpclsrc)
# acdf = pd.DataFrame(snpclsac)
# afdf = pd.DataFrame(snpclsaf)
# xticks = list(snps.iloc[:, 7].astype(str) + "_" + snps.iloc[:, 0].astype(str)+":"+snps.iloc[:, 1].astype(str) + "_" + snps['af'])
# yticks = list(rcdf.iloc[:, 0].astype(str) + "_" + rcdf.iloc[:, 1].astype(str))
# fig = plt.figure(figsize=[10, 14])
# plt.rcParams.update({'font.size': 8})
# ax = fig.add_subplot(3, 1, 1)
# ax = sns.heatmap(rcdf.iloc[:, 2:], linewidths=0.2, linecolor='k', center=0, cmap='seismic', xticklabels=False, yticklabels=yticks, cbar_kws={'fraction': 0.05})
# ax.set_title("Read coverage")
# ax = fig.add_subplot(3, 1, 2)
# ax = sns.heatmap(acdf.iloc[:, 2:], linewidths=0.2, linecolor='k', center=0, cmap='seismic', xticklabels=False, yticklabels=yticks, cbar_kws={'fraction': 0.05})
# ax.set_title("Alt allele read count")
# ax = fig.add_subplot(3, 1, 3)
# ax = sns.heatmap(afdf.iloc[:, 2:], linewidths=0.2, linecolor='k', center=0, cmap='seismic', xticklabels=xticks, yticklabels=yticks, cbar_kws={'fraction': 0.05})
# ax.set_title("Alt allele frequency")
# plt.tight_layout()
# plt.savefig(CWD+"mut_rc_iso_seq.pdf")
# plt.savefig(CWD+"mut_rc_iso_seq.png")
# plt.close()

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

    ax1 = fig.add_subplot(4, 5, i+1)
    ax1 = sns.barplot(x=list(rncnt.keys()), y=list(rncnt.values()))
    ax1.set_ylabel("Read count")
    ax1.set_ylim([0, 200000])
    ax1.tick_params(axis='x', labelrotation=90)
    ax1.ticklabel_format(axis='y', style='scientific', scilimits=(0,2))

    ax2 = fig.add_subplot(4, 5, i+2)
    ax2 = sns.barplot(x=list(rncnt.keys()), y=[rncnt[k]/genecnts[k] for k in rncnt.keys()])
    ax2.set_ylim([0, 30])
    ax2.set_ylabel("Read per gene")
    ax2.tick_params(axis='x', labelrotation=90)

    ax3 = fig.add_subplot(4, 5, i+3)
    ax3 = sns.countplot(x='clstr', data=df2)
    ax3.set_ylim([0, 200000])
    ax3.set_ylabel("Read count")
    ax3.set_xlabel(None)
    ax3.tick_params(axis='x', labelrotation=90)
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