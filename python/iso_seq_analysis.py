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

colour = {('L1', 'SNP'): '#1f77b4',
          ('L1', 'Indel'): '#84b4d6',
          ('L2', 'SNP'): '#ff7f0e',
          ('L2', 'Indel'): '#ffb97a',
          ('shared', 'SNP'): '#9467bd',
          ('shared', 'Indel'): '#c4abdb',
          ('leaf', 'SNP'): '#2ca02c',
          ('leaf', 'Indel'): '#8bcb8b',
          'L1': '#1f77b4',
          'L2': '#ff7f0e',
          'shared': '#9467bd',
          'both': '#9467bd',
          'fruit': '#c8cc92',
          'leaf': '#2ca02c',
          'theme1': '#540d6e',
          'theme2': '#ee4266',
          'theme3': '#ffd23f',
          'theme4': '#3bceac',
          'theme5': '#0ead69'}


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
    from matplotlib import pyplot as plt

    sdict = {'WT_1': 'wt1',
             'WT_19': 'wt19',
             'MUT_11_1': 'mut_11_1',
             'MUT_15': 'mut_15'}
    indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/'
    # CUTS = (10, 25, 50, 100)    # No nead to use CUTS as the new iso-seq analysis selects only good BCs
    plt.rcParams.update({'font.size': 8})
    fig=plt.figure(figsize=(6, 6))
    # Get BC read-count stats
    for i, (k, v) in enumerate(sdict.items()):
        bcrc = pd.read_table(f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/sahu_analysis/analysis/{v}_bc_readcnts.txt', header=None, delimiter=' ')
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
        # i = 1
        ax = fig.add_subplot(2, 2, i+1)
        ax = sns.regplot(x='0_x', y='0_y', data=isoillrc,
                         marker='.', color='black',
                         scatter_kws={'s': 0.5, 'color': 'lightblue'},
                         x_ci='ci', ci=95, ax=ax)

        # ax.scatter(isoillrc['0_x'], isoillrc['0_y'], s=1)
        ax.set_xlabel("scRNA read count")
        ax.set_ylabel("scIso-Seq read count")
        ax.set_title(f"Read count in {k}", fontdict={'fontsize': 8})
        # ax.set_title("Read count in cluster BC\nTotal Iso-Seq reads: {}".format(sum(isoillrc['0_y'])), fontdict={'fontsize': 8})
        # i+=1
        # ax = fig.add_subplot(3,2,i)
        # ax.hist(isoillrc['0_y'], bins=100)
        # ax.set_xlabel("scIso-Seq read count")
        # ax.set_ylabel("Number of BC")
        # i+=1
        # for c in CUTS:
        #     g = df.loc[(df[0]<1000) & (df[0]>=c)]
        #     ax = fig.add_subplot(3, 2, i)
        #     ax.set_title(i)
        #     ax.hist(g[0], bins=100)
        #     ax.set_xlim(0, 1000)
        #     ax.set_title("Number of reads: {}\nNumber of cells: {}\ncutoff: {}".format(sum(g[0]), len(g[1]), c), fontdict={'fontsize':8})
        #     i+=1
        # plt.savefig(INDIR+s+"/bc_read_dist.pdf")
        # plt.savefig(INDIR+s+"/bc_read_dist.png", transparent=True)
    plt.tight_layout(h_pad=2, w_pad=2)
    plt.savefig(indir+"/bc_read_dist.png", dpi=200, transparent=True)
    plt.close()
    return
# END
get_iso_seq_stats()



def get_allele_freq_at_sm_pos_plot():
    """
        Reads the cluster information from the scRNA analysis as done by Anshupa
        and then plots the RC and AF information for SMs in the iso-seq data
    """
    # <editor-fold desc="Define import">
    from hometools.hometools import revcomp
    from subprocess import Popen, PIPE
    import os
    from collections import defaultdict, deque
    import pandas as pd
    from itertools import product
    from matplotlib import colors as mcolors
    from matplotlib import pyplot as plt
    from scipy.cluster.hierarchy import linkage, optimal_leaf_ordering, leaves_list
    from scipy.spatial.distance import pdist
    import seaborn as sns
    # </editor-fold>


    # <editor-fold desc="Define constants and defaults">
    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/'
    sdict = {'WT_1': 'wt1',
             'WT_19': 'wt19',
             'MUT_11_1': 'mut_11_1',
             'MUT_15': 'mut_15'}
    # </editor-fold>


    # <editor-fold desc="RUN ONCE ONLY: Filter iso-seq reads from each cluster and separate them into indiviudal BAM files">
    for i, (k, v) in enumerate(sdict.items()):
        bcclt = pd.read_table(f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/sahu_analysis/analysis/{v}_bcs_clstr_id.txt', delimiter=' ')
        bcclt.columns = ['clst', 'bc']
        for grp in bcclt.groupby('clst'):
            df = grp[1]['bc'].apply(revcomp)
            df.to_csv(f'{cwd}/{k}/clstrs_{grp[0]}_bcs.txt', header=False, index=False)
            p = Popen(f'/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools view -D CB:{cwd}/{k}/clstrs_{grp[0]}_bcs.txt -O BAM -o {cwd}/{k}/clstrs_{grp[0]}_reads.bam {cwd}/{k}/{k}.CB_sorted.bam'.split(), stdout=PIPE, stderr=PIPE)
            o = p.communicate()
            if o[1] != b'':
                print(f"error in get reads for sample: {k}, cluster: {grp[0]}")
            p = Popen(f'/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools fasta -0 {cwd}/{k}/clstrs_{grp[0]}_reads.fa.gz -T CB {cwd}/{k}/clstrs_{grp[0]}_reads.bam'.split(), stdout=PIPE, stderr=PIPE)
            o = p.communicate()
    # </editor-fold>

    # Align the above reads and get read-counts at SM positions
    # using code in iso_seq_analysis.sh

    return
# END


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