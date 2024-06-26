import argparse


# Set colours
colour = {('L1', 'SNV'): '#1f77b4',
          ('L1', 'Indel'): '#84b4d6',
          ('L2', 'SNV'): '#ff7f0e',
          ('L2', 'Indel'): '#ffb97a',
          ('shared', 'SNV'): '#9467bd',
          ('shared', 'Indel'): '#c4abdb',
          ('fruit', 'SNV'): '#003399',
          ('fruit', 'Indel'): '#6699ff',
          ('leaf', 'SNV'): '#2ca02c',
          ('leaf', 'Indel'): '#8bcb8b',
          'L1': '#1f77b4',
          'L2': '#ff7f0e',
          'L1+leaf': '#34C39C',
          'L2+leaf': '#D1C929',
          'shared': '#9467bd',
          'both': '#9467bd',
          'fruit': '#003399',
          'leaf': '#2ca02c',
          'any': '#540d6e',
          'theme1': '#540d6e',
          'theme2': '#ee4266',
          'theme3': '#ffd23f',
          'theme4': '#3bceac',
          'theme5': '#0ead69'}


def read_depth_summary(rd, chrsize, output, min=50):
    from matplotlib import pyplot as plt
    import pandas as pd
    import numpy as np
    contigs_depth_sum = {}
    with open(rd, 'r') as f:
        contig = ''
        sum = 0
        FIRST = True
        # count = 0
        for line in f:
            line = line.strip().split()
            if line[0] != contig:
                if not FIRST: contigs_depth_sum[contig] = sum
                else: FIRST = False
                contig  = line[0]
                sum     = int(line[2])
                print(contig)
            else:
                sum += int(line[2])
        contigs_depth_sum[contig] = sum

    df = pd.read_table(chrsize, header=None)
    contigs_sizes = dict(zip(df[0], df[1]))
    contigs_list = list(contigs_sizes.keys())

    contigs_mean_depth = {k:[contigs_depth_sum[k]/contigs_sizes[k] if k in contigs_depth_sum else 0][0] for k in contigs_list}
    df['md'] = [contigs_mean_depth[k] for k in df[0]]
    df.sort_values(['md'], ascending=True, inplace=True)
    df['cs'] = np.cumsum(df[1])

    fig = plt.figure(figsize=[6,10])
    ax1 = fig.add_subplot(3,1,1)
    ax1.scatter(df.md, df[1])
    ax1.set_xlim([0,min])
    ax1.set_xlabel("mean read depth")
    ax1.set_ylabel("contig size")

    ax2 = fig.add_subplot(3,1,2)
    plt.hist(df.md, bins=min, range=[0,min])
    ax2.set_xlabel("mean read depth")
    ax2.set_ylabel("Number of contigs")
    ax2.set_xlim([0,min])

    ax3 = fig.add_subplot(3,1,3)
    plt.plot(df.md, df.cs)
    ax3.set_xlim([0, min])
    ax3.set_xlabel("mean read depth")
    ax3.set_ylabel("Cumulative assembly size")
    ax3.minorticks_on()
    ax3.grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.25)
    ax3.grid(b=True, which='major', axis='both', linestyle='-', linewidth=2)
    plt.subplots_adjust(left=0.1, bottom=0.05, right=0.97, top=0.98, wspace=0.1, hspace=0.25)
    fig.savefig(output)
    return df


def read_coverage(depthfin, fout, mm2_good, y_cut=10):
    from collections import deque
    import numpy as np
    from matplotlib import pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.patches import Rectangle

    cnts = deque()
    chr_dep = {}

    # depthfin = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/candidate.contigs.v3.unmapped_cursr.sorted.depth'
    with open(depthfin, 'r') as f:
        chr = ''
        for line in f:
            c, _, d = line.strip().split()
            if c == chr:
                cnts.append(int(d))
            elif chr == '':
                chr = c
                cnts.append(int(d))
            else:
                chr_dep[chr] = {}
                cnts = np.array(cnts)
                for i in range(0, len(cnts), 10000):
                    chr_dep[chr][(i + min(i+50000, len(cnts)))/2] = np.mean(cnts[i: min(i+50000, len(cnts))])
                cnts = deque()
                chr = c
                cnts.append(int(d))
        chr_dep[chr] = {}
        cnts = np.array(cnts)
        for i in range(0, len(cnts), 10000):
            chr_dep[chr][(i + min(i+50000, len(cnts)))/2] = np.mean(cnts[i: min(i+50000, len(cnts))])



    # mm2_good = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/minimap2_check/contigs.low_overlap.txt'
    mm2_out = {}
    with open(mm2_good, 'r') as f:
        for line in f:
            line = line.strip().split()
            try:
                mm2_out[line[0]] = [float(line[1]), line[2]]
            except IndexError as e:
                mm2_out[line[0]] = [float(line[1]), ""]
    # mm2_ks = list(mm2_out.keys())
    mm2_ks = sorted(mm2_out.keys(), key=lambda x: max(chr_dep[x].keys()), reverse=True)

    # with PdfPages('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/TMP_contig_depths.pdf') as pdf:
    with PdfPages(fout) as pdf:
        count = 0
        fig = plt.figure(figsize=[8, 16])
        index = 0
        for k in mm2_ks:
            count += 1
            index += 1
            ax = fig.add_subplot(10, 1, index)
            ax.plot(list(chr_dep[k].keys()), list(chr_dep[k].values()), c='black')
            ax.axhline(y=y_cut, color='grey', linestyle='--')
            if mm2_out[k][0] < 0.5:
                ax.spines['top'].set_color("red")
                ax.spines['top'].set_linewidth(3)
                ax.spines['right'].set_color("red")
                ax.spines['right'].set_linewidth(3)
                ax.spines['bottom'].set_color("red")
                ax.spines['bottom'].set_linewidth(3)
                ax.spines['left'].set_color("red")
                ax.spines['left'].set_linewidth(3)
            ax.set_ylabel(k)
            plt.ylim(0, max(50, max(list(chr_dep[k].values()))))

            for j in mm2_out[k][1].split(';'):
                try:
                    s, e = j.split('-')
                except ValueError as e:
                    continue
                s, e = int(s), int(e)

                ax.add_patch(Rectangle((s, 0), e-s, max(50, max(list(chr_dep[k].values()))),
                                       edgecolor='lightgrey',
                                       facecolor='lightgrey',
                                       fill=True,
                                       linewidth=0.1))



            if count%10 == 0:
                index = 0
                plt.subplots_adjust(left=0.1, bottom=0.05, right=0.9, top=0.98, wspace=0.1, hspace=0.2)
                pdf.savefig()
                plt.close()
                fig = plt.figure(figsize=[8, 16])
                # break
        plt.subplots_adjust(left=0.1, bottom=0.05, right=0.9, top=0.98, wspace=0.1, hspace=0.2)
        pdf.savefig()
        plt.close()


def assembly_plots():
    """
        Collection of functions for anlysing and visualising assembly quality
    """
    # <editor-fold desc="Define imports">
    from itertools import permutations
    from collections import deque, Counter, defaultdict, OrderedDict
    import subprocess
    import numpy as np
    import pandas as pd
    from matplotlib import pyplot as plt
    import seaborn as sns
    import igraph as ig
    import pybedtools as bt
    from scipy.cluster.hierarchy import dendrogram, linkage, optimal_leaf_ordering, leaves_list
    from scipy.spatial.distance import squareform
    from scipy.sparse.csgraph import minimum_spanning_tree
    from scipy.sparse import triu, csr_matrix
    from hometools.plot import plot2emf, cleanax
    from hometools.hometools import printdf, unlist, canonical_kmers, Namespace, cntkmer, revcomp, readfasta, mergeRanges, sumranges
    from multiprocessing import Pool
    import json
    from itertools import product
    # </editor-fold>

    # <editor-fold desc="Define default values">
    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemblyplots/'
    curfa = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
    orafa = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/ora.genome.v1.fasta'
    curlen = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.chrlen.txt'
    oralen = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/ora.genome.chrlen.txt'
    cur5tel = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/telomer.5end.kmer.cur.v1.bed'
    cur3tel = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/telomer.3end.kmer.cur.v1.bed'
    ora5tel = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/telomer.5end.kmer.ora.v1.bed'
    ora3tel = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/telomer.3end.kmer.ora.v1.bed'
    curgff = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/EVM_PASA/pasa_on_mancur/cur.pasa_out.sort.protein_coding.3utr.gff3'
    oragff = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/ora/EVM_PASA/pasa_on_mancur/ora.pasa_out.sort.protein_coding.3utr.gff3'

    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    # </editor-fold>

    # <editor-fold desc="Plot positions of telomeric repeats on the chromosomes>
    def gethist(cl, ct, maxx, n, t):
        fig = plt.figure(figsize= (4, 7), dpi=300)
        for i in range(0, 8):
            ax = fig.add_subplot(8, 1, i+1)
            ax.spines[:].set_visible(False)
            ax.set_xlim([0, maxx])
            c = cl.iloc[i, 0]
            l = cl.iloc[i, 1]
            ax.hlines(0, 0, l, linewidth=5)
            bins = np.arange(0, l, 200000)
            h = np.histogram(ct.loc[ct[0] == c, 1], bins=bins)
            ys = h[0]
            xs = bins[1:] - 100000
            ax.bar(x=xs,height=ys, width=180000)
            ax.set_xlabel(c)
        fig.supylabel(f'count')
        fig.suptitle(t)
        plt.tight_layout()
        plt.savefig(n)
        plt.close()
        return
    # END

    # get plot for currot
    clen = pd.read_table(curlen, header=None).iloc[:8]
    maxx = max(clen[1])
    # 5' telomers
    tel5 = pd.read_table(cur5tel, header=None)
    gethist(clen, tel5, maxx, f'{cwd}/cur.telo.5end.png', "Currot 5' telo (3x) sequence")
    # 3' telomers
    tel3 = pd.read_table(cur3tel, header=None)
    gethist(clen, tel3, maxx, f'{cwd}/cur.telo.3end.png', "Currot 3' telo (3x) sequence")

    # get plot for orangered
    olen = pd.read_table(oralen, header=None).iloc[:8]
    maxx = max(olen[1])
    # 5' telomers
    tel5 = pd.read_table(ora5tel, header=None)
    gethist(olen, tel5, maxx, f'{cwd}/ora.telo.5end.png', "Orangered 5' telo (3x) sequence")
    # 3' telomers
    tel3 = pd.read_table(ora3tel, header=None)
    gethist(olen, tel3, maxx, f'{cwd}/ora.telo.3end.png', "Orangered 3' telo (3x) sequence")

    # </editor-fold>


    # <editor-fold desc="Plot gene distribution on the chromosomes">
    # Get currot gene distribution
    clen = pd.read_table(curlen, header=None).iloc[:8]
    maxx = max(clen[1])
    cgff = pd.read_table(curgff, header=None)
    cgff = cgff.loc[cgff[2] == 'gene', [0, 3]]
    cgff.columns = [0, 1]
    gethist(clen, cgff, maxx, f'{cwd}/currot.gene.distribution.png', "Currot gene distribution")
    # Get gene positions
    olen = pd.read_table(oralen, header=None).iloc[:8]
    maxx = max(olen[1])
    ogff = pd.read_table(oragff, header=None)
    ogff = ogff.loc[ogff[2] == 'gene', [0, 3]]
    ogff.columns = [0, 1]
    gethist(olen, ogff, maxx, f'{cwd}/orangered.gene.distribution.png', "Orangered gene distribution")

    # </editor-fold>
    return
# END


def mutation_spectra():
    """
    Plots describing mutations spectra of somatic mutations identified from the
    layer_specific sequencing
    """
    # Figure size: Paper: The width of figures, when printed, will usually be 5.7 cm (2.24 inches or 1 column), 12.1 cm (4.76 inches or 2 columns), or 18.4 cm (7.24 inches or 3 columns).

    # <editor-fold desc="Define imports">
    from itertools import permutations
    from collections import deque, Counter, defaultdict, OrderedDict
    from subprocess import Popen, PIPE
    import numpy as np
    import pandas as pd
    from matplotlib import pyplot as plt
    from matplotlib import cm as mcm
    import matplotlib.transforms as mtransforms
    from matplotlib_venn import venn3, venn2
    import matplotlib.image as mpimg
    import seaborn as sns
    import igraph as ig
    import pybedtools as bt
    from scipy.cluster.hierarchy import dendrogram, linkage, optimal_leaf_ordering, leaves_list, set_link_color_palette
    from scipy.spatial.distance import squareform, pdist
    from scipy.sparse.csgraph import minimum_spanning_tree
    from scipy.sparse import triu, csr_matrix
    from scipy.stats import ttest_rel, fisher_exact, spearmanr
    from hometools.plot import plot2emf, cleanax, inlineplotsr
    from hometools.hometools import printdf, unlist, canonical_kmers, Namespace, cntkmer, revcomp, readfasta, mergeRanges, sumranges, readsyriout, getvalues
    from multiprocessing import Pool
    import json
    from itertools import product
    from statsmodels.stats.multitest import multipletests
    from statsmodels.stats.weightstats import DescrStatsW
    from tqdm import tqdm
    import os
    import pyranges as pr
    import re
    from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
    from Bio.Phylo import draw as pdraw
    from Bio.Phylo import to_networkx
    import networkx
    from itertools import combinations
    # </editor-fold>


    # <editor-fold desc="Define default values">
    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
    refcur = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
    curfai = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.fai'
    syriout = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out'
    syrisnppath = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.snps.txt'
    syriindpath = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.indels.txt'
    allsm = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_in_all_samples.manually_selected.cleaned.noisy_annotated.csv'
    allsmrcdir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_in_all_samples_readcounts/'
    allsmrc = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_in_all_samples.read_counts.txt'
    smbcs = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/scrna_clusters/bcs_with_sm_reads.txt'
    figout = '/local/goel/Dropbox/projects/apricot_leaf/manuscript/figures/'
    scrnadir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/scrna_clusters/'
    scisodir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/'
    branches = ['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']
    branchorder = ['wt_7', 'wt_1', 'wt_18', 'wt_19', 'mut_15', 'mut_11_1', 'mut_11_2']
    bnamedict = dict(zip(branchorder, [f'B{i}' for i in range(1, 8)]))
    SMALL_SIZE = 7
    BIGGER_SIZE = 8
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=SMALL_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=SMALL_SIZE)  # fontsize of the figure title
    plt.rc('font', family='Arial')  # fontsize of the figure title
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42
    scovs = {'wt_1':  {'l1': (40, 240), 'l2': (20, 180)},
             'wt_7':  {'l1': (30, 180), 'l2': (30, 180)},
             'wt_18': {'l1': (40, 200), 'l2': (30, 180)},
             'wt_19': {'l1': (40, 220), 'l2': (40, 220)},
             'mut_11_1': {'l1': (20, 160), 'l2': (20, 150)},
             'mut_11_2': {'l1': (30, 200), 'l2': (20, 180)},
             'mut_15': {'l1': (40, 220), 'l2': (30, 220)}}
    genomesize = 233151598
    chrom_genome_size = 228047739
    # </editor-fold>


    # <editor-fold desc="Read the layer_specific SM as well as SMs shared between layers and generate one layer SM object">
    layersp = pd.read_table(f'{cwd}all_layer_specific_somatic_variants.txt')    # layer_specific_variations
    ## Because reads L3 are not enriched for L3 and are mostly L2, removing them
    layersp.loc[layersp.Layer == 'L1,L3', 'Layer'] = 'L1'
    layersp.loc[layersp.Layer == 'L2,L3', 'Layer'] = 'L2'
    nrc = deque()
    naf = deque()
    for row in layersp.itertuples(index=False):
        if np.isnan(row.read_count):
            l = row.Layer
            rc = row[row._fields.index(f'{l}_RC')]
            af = row[row._fields.index(f'{l}_AF')]
            nrc.append(rc)
            naf.append(af)
        else:
            nrc.append(row.read_count)
            naf.append(row.allele_freq)
    layersp.read_count = nrc
    layersp.allele_freq = naf

    sharedsm = pd.read_table(f'{cwd}/all_samples_candidate_high_cov_sms3.only_selected.csv', header=None)
    sharedsm.columns = ['chromosome', 'position', 'branch', 'Layer', 'ref_allele', 'alt_allele', 'read_count', 'allele_freq', 'Selected', 'Ploidy', 'Remarks']
    sharedsm.Layer = sharedsm.Layer.str.upper()
    # TODO: ADD THESE POSITIONS TO GENE CONVERSIONS
    # filtering positions that are present in all but one sample as it is not clear whether they are de novo mutation or gene conversions
    sharedsm = sharedsm.loc[sharedsm.position != 367749]
    sharedsm = sharedsm.loc[sharedsm.position != 16736100]

    data = pd.concat([layersp, sharedsm])
    data['type'] = 'SNP'
    data.loc[[a[0] in '+-' for a in data.alt_allele], 'type'] = 'Indel'
    datafilter = data.loc[:, ['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele', 'read_count', 'allele_freq', 'Layer', 'type']].drop_duplicates(subset=['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele', 'Layer', 'type']).copy()    # Filtered dataframe to be used for generating statistics
    datafilter.sort_values('chromosome position'.split(), inplace=True)
    # datafilter.to_csv('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.filtered.txt', index=False, header=True, sep='\t')

    # </editor-fold>


    # <editor-fold desc="OBSOLETE: Read and filter layer SM calls">

    datafilter = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.filtered.txt')
    # Annotate SMs as L1, L2 or shared
    l1sms = datafilter.loc[datafilter.Layer == 'L1', 'chromosome position ref_allele alt_allele'.split()].to_numpy()
    l1pos = set(list(map(tuple, l1sms)))
    l2sms = datafilter.loc[datafilter.Layer == 'L2', 'chromosome position ref_allele alt_allele'.split()].to_numpy()
    l2pos = set(list(map(tuple, l2sms)))
    sharedpos = l1pos.intersection(l2pos)
    smtype = ['shared' if (row[0], row[1], row[3], row[4]) in sharedpos else row[7] for row in datafilter.itertuples(index=False)]
    datafilter.Layer = smtype
    datafilter.drop_duplicates(subset=['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele', 'Layer', 'type'], inplace=True)
    # </editor-fold>


    # <editor-fold desc="Read and filter all somatic mutation file">
    allsmdata = pd.read_csv(allsm)
    datafilter = allsmdata.loc[allsmdata.tissue.isin('L1 L2'.split())]
    datafilter.columns = 'chromosome position alt_allele branch Layer organ Low_confidence Potential_alt_allele'.split()
    refseq = readfasta(refcur)
    smpos = set(map(tuple, datafilter[['chromosome', 'position']].to_numpy()))
    ref_allele = deque()
    for p in smpos:
        ref_allele.append([p[0], p[1], refseq[p[0]][p[1]-1: p[1]]])
    ref_allele = pd.DataFrame(ref_allele)
    ref_allele.columns = 'chromosome position ref_allele'.split()
    datafilter = datafilter.merge(ref_allele)
    datafilter = datafilter.loc[:, 'chromosome position branch ref_allele alt_allele Layer'.split()]
    # Annotate SMs as L1, L2 or shared
    l1sms = datafilter.loc[datafilter.Layer == 'L1', 'chromosome position ref_allele alt_allele'.split()].to_numpy()
    l1pos = set(list(map(tuple, l1sms)))
    l2sms = datafilter.loc[datafilter.Layer == 'L2', 'chromosome position ref_allele alt_allele'.split()].to_numpy()
    l2pos = set(list(map(tuple, l2sms)))
    sharedpos = l1pos.intersection(l2pos)
    smtype = ['shared' if (row[0], row[1], row[3], row[4]) in sharedpos else row[5] for row in datafilter.itertuples(index=False)]
    datafilter.Layer = smtype
    datafilter['type'] = 'SNV'
    datafilter.loc[[a[0] in '+-' for a in datafilter.alt_allele], 'type'] = 'Indel'
    datafilter.drop_duplicates(subset=['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele', 'Layer', 'type'], inplace=True)
    # </editor-fold>


    # <editor-fold desc="Get stats/counts">

    layerjsondata = dict()
    layerjsondata['Number of SM events'] = datafilter.shape[0]
    layerjsondata['Number of unique SNPs'] = datafilter.loc[datafilter.type == "SNP"].drop_duplicates(subset=["chromosome", "position", "ref_allele", "alt_allele"]).shape[0]
    layerjsondata['Number of unique Indel'] = datafilter.loc[datafilter.type == "Indel"] \
        .drop_duplicates(subset=["chromosome", "position", "ref_allele", "alt_allele"]) \
        .shape[0]
    layerjsondata['Total SM counts'] = datafilter \
        .drop_duplicates(subset=["chromosome", "position", "alt_allele"]) \
        .shape[0]
    layerjsondata['Number of SMs in L1'] = datafilter.loc[datafilter.Layer == "L1"].drop_duplicates(subset=["chromosome", "position", "ref_allele", "alt_allele"]).shape[0]
    layerjsondata['Number of SMs in L2'] = datafilter.loc[datafilter.Layer == "L2"].drop_duplicates(subset=["chromosome", "position", "ref_allele", "alt_allele"]).shape[0]
    layerjsondata['Number of shared SMs'] = datafilter.loc[datafilter.Layer == "shared"].drop_duplicates(subset=["chromosome", "position", "ref_allele", "alt_allele"]).shape[0]
    layerjsondata['Number of branches for each SM'] = OrderedDict(sorted(Counter([len(set(grp[1].branch)) for grp in datafilter.groupby(['chromosome', 'position', 'alt_allele'])]).items()))
    layerjsondata["Number of SMs per branch"] = {grp[0]: grp[1].drop_duplicates(subset=["chromosome", "position", "ref_allele", "alt_allele"]).shape[0] for grp in datafilter.groupby('branch')}
    # Is the statistically significance? No
    g = datafilter.branch.value_counts()
    z = (g - np.mean(g))/np.std(g)  # Z-transform shows that No value is outside the 95% range
    layerjsondata["Dimeric indels"] = {i: Counter([grp[0][2] for grp in datafilter.loc[datafilter.type == 'Indel'].groupby('chromosome position alt_allele'.split())])[i] for i in ['-AT', '-TA', '+AT', '+TA']}
    layerjsondata["Dimeric indels total count"] = sum([Counter([grp[0][2] for grp in datafilter.loc[datafilter.type == 'Indel'].groupby('chromosome position alt_allele'.split())])[i] for i in ['-AT', '-TA', '+AT', '+TA']])

    layerjsondata["Number of L1 SMs in all branches"] = Counter(datafilter.loc[datafilter.Layer == 'L1'].groupby('chromosome position alt_allele'.split()).alt_allele.value_counts().to_list())[7]
    layerjsondata["Number of L2 SMs in all branches"] = Counter(datafilter.loc[datafilter.Layer == 'L2'].groupby('chromosome position alt_allele'.split()).alt_allele.value_counts().to_list())[7]
    layerjsondata["Number of shared SMs in all branches"] = Counter(datafilter.loc[datafilter.Layer == 'shared'].groupby('chromosome position alt_allele'.split()).alt_allele.value_counts().to_list())[7]
    with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/axillary_plots/layer_sm_stats.json', 'w') as fout:
        json.dump(layerjsondata, fout, indent=2)

    # </editor-fold>

    def pretty_ticks(xmin, xmax, nticks, ax, axis='x'):
        ticks = np.linspace(xmin, xmax, nticks)
        if axis == 'x':
            ax.set_xticks(ticks)
        elif axis == 'y':
            ax.set_yticks(ticks)
        return ax


    # <editor-fold desc="Figure 1 plots">
    barwidth = 0.125
    dpi = 100
    pltwidth = 7.24
    fig1 = plt.figure(figsize=[pltwidth, 8.5], dpi=dpi)
    trans = mtransforms.ScaledTranslation(-20/dpi, 5/dpi, fig1.dpi_scale_trans)


    # <editor-fold desc="A: Add tree image">
    def getcntplt(ax, branch, bar_width):
        # ax = fig1.add_axes(axispos)
        y_bottom = [0]
        ymax = 75
        ax.set_xlim([-0.5, 0.5])
        ax.set_ylim([0, ymax])
        ax.set_yticks([0, ymax])
        ax.tick_params(bottom=False, labelbottom=False)
        ax.set_xlabel(bnamedict[branch], fontdict={'backgroundcolor': 'white', 'linespacing': 0, 'horizontalalignment': 'center', 'bbox': {'pad': 1, 'fc': 'white', 'ec': 'white'}})
        ax.set_yticklabels(labels=ax.get_yticklabels(), fontdict={'backgroundcolor': 'white', 'linespacing': 0, 'horizontalalignment': 'center', 'verticalalignment': 'center', 'bbox': {'pad': 1, 'fc': 'white', 'ec': 'white'}})
        ax.spines[:].set_visible(False)
        for l in ('L1', 'L2'):
            for t in ('SNV', 'Indel'):
                y = branchscnt[(l, t)][branch]
                ax.bar([0], [y], label=f'{l} {t}', bottom=y_bottom, color=colour[(l, t)], width=bar_width)
                y_bottom = [y_bottom[0] + y]
        return ax
    # END

    tree = mpimg.imread(f'{cwd}/../tree_sketch_fruit_nolength.png')
    axA = fig1.add_axes([0.03, 0.675, 0.4, 0.3])
    axA.imshow(tree, aspect='auto')
    axA.text(0, 1, 'A', transform=axA.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    axA.spines[:].set_visible(False)
    axA.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)

    branchscnt = {g[0]: defaultdict(int, g[1].branch.value_counts()) for g in datafilter.groupby(['Layer', 'type'])}
    bar_width = 0.75
    axwidth = (barwidth*1/pltwidth)/0.4
    axA1 = getcntplt(axA.inset_axes([0.9, 0.7, axwidth, 0.1]), 'wt_1', bar_width)
    axA2 = getcntplt(axA.inset_axes([0.9, 0.25, axwidth, 0.1]), 'wt_7', bar_width)
    axA3 = getcntplt(axA.inset_axes([0.75, 0.85, axwidth, 0.1]), 'wt_18', bar_width)
    axA4 = getcntplt(axA.inset_axes([0.45, 0.85, axwidth, 0.1]), 'wt_19', bar_width)
    axA5 = getcntplt(axA.inset_axes([0.25, 0.85, axwidth, 0.1]), 'mut_15', bar_width)
    axA6 = getcntplt(axA.inset_axes([0.1, 0.85, axwidth, 0.1]), 'mut_11_1', bar_width)
    axA7 = getcntplt(axA.inset_axes([0.1, 0.45, axwidth, 0.1]), 'mut_11_2', bar_width)

    # </editor-fold>


    # <editor-fold desc="B: Number of SMs per layer">
    axwidth = barwidth*3/pltwidth
    axB = fig1.add_axes([0.5, 0.875, axwidth, 0.1])
    y_bottom = [0]*3
    # y_bottom = [0]*2
    for t in ('SNV', 'Indel'):
        layercnts = datafilter.loc[datafilter.type == t] \
            .drop_duplicates(subset='chromosome position alt_allele'.split()) \
            .Layer.value_counts().to_dict()
        axB.bar(['L1', 'L2', 'shared'], [layercnts[l] for l in ('L1', 'L2', 'shared')], bottom=y_bottom, color=[colour[l, t] for l in ('L1', 'L2', 'shared')], width=0.75)
        y_bottom = [layercnts[l] for l in ('L1', 'L2', 'shared')]
    axB = pretty_ticks(0, 125, 6, axB, 'y')
    axB.set_xlabel('')
    axB.set_ylabel("SMs")
    axB.tick_params(axis='x', rotation=45)
    axB = cleanax(axB)
    axB.text(0, 1, 'B', transform=axB.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')

    # </editor-fold>


    # # <editor-fold desc="C: Number of branches covered by each SM">
    #
    # pos = datafilter.drop_duplicates(subset=['chromosome', 'position', 'alt_allele']).copy()
    # samplecnt = {grp[0]: len(set(grp[1].branch)) for grp in datafilter.groupby(['chromosome', 'position', 'alt_allele'])}
    # pos['Branch Count'] = [samplecnt[row.chromosome, row.position, row.alt_allele] for row in pos.itertuples(index=False)]
    # n_branch = pos.groupby(['Layer', 'type', 'Branch Count']).size().to_dict()
    # axwidth = barwidth*7/pltwidth
    # axC = fig1.add_axes([0.625, 0.875, axwidth, 0.1])
    # y_bottom = [0]*7
    # for l in ('L1', 'L2'):
    #     for t in ('SNV', 'Indel'):
    #         y = [n_branch[(l, t, i)] if (l, t, i) in n_branch else 0 for i in range(1, 8)]
    #         axC.bar(np.arange(0, 7), y, label=f'{l}_{t}', bottom=y_bottom, color=colour[(l, t)], width=0.75)
    #         y_bottom = [y_bottom[i] + y[i] for i in range(7)]
    # axC.set_xlim(-0.5, 6.5)
    # axC.set_xlabel("Number of branch")
    # axC.set_ylabel("SMs")
    # axC.set_xticks(ticks=range(7))
    # axC.set_xticklabels(np.arange(1, 8))
    # axC = cleanax(axC)
    # axC.text(0, 1, 'C', transform=axC.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    #
    # # </editor-fold>


    # <editor-fold desc="F: Plot mutation load">
    """Consider each sample as separate lineage and divide by the diploid genome size"""
    # mrates = allsmdata.groupby('branch tissue'.split()).position.nunique()/480000000
    mrates = allsmdata.groupby('branch tissue'.split()).position.nunique()
    for b, l in mrates.index.values:
        if l == 'leaf':
            continue
        uncallable = bt.BedTool(f'{bnamedict[b]}_{l.lower()}_uncallable_positions.merged.bed').filter(lambda x: 'CUR' in x[0]).saveas()
        uncallable_size = sum([int(b[2]) - int(b[1]) for b in uncallable])
        mrates[(b, l)] /= (genomesize - uncallable_size)*2
    mrates = pd.DataFrame(mrates)
    mrates.reset_index(inplace=True)
    mrates.columns = 'branch tissue value'.split()
    mrates['mrate'] = 'Mutation load'
    mrates = mrates.loc[mrates.tissue !='leaf']
    axF = fig1.add_axes([0.635, 0.875, 0.075,  0.1])
    axF = sns.violinplot(data=mrates, x='mrate', y="value", color=colour['theme1'], linewidth=0, inner=None, ax=axF)
    plt.setp(axF.collections, alpha=.2)
    ylim = axF.get_ylim()
    axF = sns.stripplot(data=mrates, x="mrate", y="value", jitter=True, palette=colour, zorder=1, ax=axF, hue='tissue', legend=False, s=3)
    axF = pretty_ticks(0, 2e-07, 3, axF, 'y')
    axF.set_ylim(ylim[0], 2e-07)
    axF.set_title('')
    axF.set_ylabel('Mutation load (x$10^{-7}$)')
    axF.set_xlabel("All SMs")
    axF.ticklabel_format(axis='y', useOffset=False, style='plain')
    yticks = axF.get_yticks()
    yticksl = yticks*10000000
    axF.set_yticklabels(yticksl)
    # axF.set_yticks(yticks[1:-1])
    # axF.set_yticklabels(yticksl[1:-1])
    axF.tick_params(axis='x', bottom=False, labelbottom=False)
    axF = cleanax(axF)
    axF.text(0, 1, 'C', transform=axF.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')

    # mutsinallL1 = Counter(allsmdata.loc[allsmdata.tissue == 'L1'].groupby('position alt_allele'.split()).position.count())[7]       # Number of SMs that are in all branches
    # mrates = allsmdata.groupby('branch tissue'.split()).position.nunique()
    # mrates = pd.DataFrame(mrates)
    # mrates.reset_index(inplace=True)
    # mrates.columns = 'branch tissue value'.split()
    # mrates.loc[mrates.tissue == 'L1', 'value'] -= mutsinallL1
    # mrates['value'] /= 480000000
    # mrates['mrate'] = 'Mutation load'

    # axE2 = fig3.add_axes([0.225, 0.075, 0.1, 0.3])
    # axE2 = sns.violinplot(data=mrates, x='mrate', y="value", color=colour['theme1'], linewidth=0, inner=None, ax=axE2)
    # plt.setp(axE2.collections, alpha=.2)
    # ylim = axE2.get_ylim()
    # axE2 = sns.stripplot(data=mrates, x="mrate", y="value", jitter=True, palette=colour, zorder=1, ax=axE2, hue='tissue', legend=False, size=3)
    # axE2.set_ylim(ylim)
    # axE2.set_ylabel('')
    # # axE2.set_xlabel("Filtered SMs present\nin L1 of branches")
    # axE2.set_xlabel("Filtered SMs")
    # axE2.ticklabel_format(axis='y', useOffset=False, style='plain')
    # axE2.set_ylim(axE1.get_ylim())
    # plt.tick_params(axis='x', bottom=False, labelbottom=False)
    # plt.tick_params(axis='y', left=False, labelleft=False)
    # axE2 = cleanax(axE2)
    # axE2.spines['left'].set_visible(False)
    # axE.append(axE2)


    # </editor-fold


    # # <editor-fold desc="D: Plot for number of SMs per branch">
    #
    # branchscnt = {g[0]: defaultdict(int, g[1].branch.value_counts()) for g in datafilter.groupby(['Layer', 'type'])}
    # axD = fig1.add_axes([0.55, 0.65, 0.45, 0.1])
    # y_bottom = [0]*7
    # bar_width = 0.5
    # # for l in ('L1', 'L2', 'shared'):
    # for l in ('L1', 'L2'):
    #     for t in ('SNP', 'Indel'):
    #         y = [branchscnt[(l, t)][b] for b in branches]
    #         axD.bar(np.arange(1, 8), y, label=f'{l} {t}', bottom=y_bottom, color=colour[(l, t)], width=bar_width)
    #         y_bottom = [y_bottom[i] + y[i] for i in range(7)]
    # axD.set_xlabel("Branch")
    # axD.set_ylabel("Number of SMs")
    # axD.set_xticks(np.arange(1, 8))
    # axD.set_xticklabels(branches)
    # axD.spines['right'].set_visible(False)
    # axD.spines['top'].set_visible(False)
    # # axD.legend(bbox_to_anchor=(1.01, 1), frameon=False)
    # # axD.tick_params(axis='x', rotation=45)
    # axD.text(0, 1, 'D', transform=axD.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    #
    branchscnt = {g[0]: defaultdict(int, g[1].branch.value_counts()) for g in datafilter.groupby('Layer')}
    garb1 = [branchscnt['L1'][b] - 21 for b in branches]
    garb2 = [branchscnt['L2'][b] - 1 for b in branches]
    ttest_rel(garb1, garb2, alternative='greater')
    # # </editor-fold>


    # <editor-fold desc="C: Get transversion vs transitions plot">
    unisnp = datafilter.loc[datafilter.type == 'SNV'].drop_duplicates(subset=['chromosome', 'position', 'ref_allele', 'alt_allele', 'Layer'])
    unisnp['snptype'] = 'transversions'
    unisnp.loc[(unisnp['ref_allele'] == 'A') & (unisnp['alt_allele'] == 'G'), 'snptype'] = 'transitions'
    unisnp.loc[(unisnp['ref_allele'] == 'C') & (unisnp['alt_allele'] == 'T'), 'snptype'] = 'transitions'
    unisnp.loc[(unisnp['ref_allele'] == 'G') & (unisnp['alt_allele'] == 'A'), 'snptype'] = 'transitions'
    unisnp.loc[(unisnp['ref_allele'] == 'T') & (unisnp['alt_allele'] == 'C'), 'snptype'] = 'transitions'
    mutchange = list(permutations('ATGC', 2))
    mutchange_counts = dict()
    for grp in unisnp.groupby('Layer'):
        mutchange_counts[grp[0]] = dict()
        for m in mutchange:
            mutchange_counts[grp[0]]['>'.join(m)] = grp[1].loc[(grp[1].ref_allele == m[0]) & (grp[1].alt_allele == m[1])].shape[0]
    df = pd.DataFrame(mutchange_counts)
    df['var'] = df.index.values
    df.reset_index(drop=True, inplace=True)
    df.loc[:, 'L1'] = round(df.loc[:, 'L1']*100/sum(df.L1), 2)
    df.loc[:, 'L2'] = round(df.loc[:, 'L2']*100/sum(df.L2), 2)
    df.loc[:, 'shared'] = round(df.loc[:, 'shared']*100/sum(df.shared), 2)
    df = df.melt(id_vars=['var'])
    df.columns = ['SNP', 'Layer', 'p']
    garb = deque()
    for grp in df.groupby('Layer'):
        d = grp[1]
        # print(d)
        garb.append(['A,T>G,C', grp[0], sum(d.loc[d.SNP.isin(('A>G', 'T>C')), 'p'])])
        garb.append(['C,G>T,A', grp[0], sum(d.loc[d.SNP.isin(('G>A', 'C>T')), 'p'])])
        garb.append(['A,T>C,G', grp[0], sum(d.loc[d.SNP.isin(('A>C', 'T>G')), 'p'])])
        garb.append(['A,T>T,A', grp[0], sum(d.loc[d.SNP.isin(('A>T', 'T>A')), 'p'])])
        garb.append(['C,G>A,T', grp[0], sum(d.loc[d.SNP.isin(('C>A', 'G>T')), 'p'])])
        garb.append(['C,G>G,C', grp[0], sum(d.loc[d.SNP.isin(('C>G', 'G>C')), 'p'])])
    df = pd.DataFrame(garb)
    df.columns = ['SNP type', 'Layer', 'Percentage of SMs']
    df = df.loc[~(df.Layer == 'shared')]
    axwidth = barwidth*12/pltwidth
    # axB = fig1.add_axes([0.5, 0.875, axwidth, 0.1])
    axE = fig1.add_axes([0.78, 0.875, axwidth, 0.1])
    axE = sns.barplot(data=df, x='SNP type', y='Percentage of SMs', hue='Layer', palette=colour, order=['A,T>G,C', 'C,G>T,A', 'A,T>C,G', 'A,T>T,A', 'C,G>A,T', 'C,G>G,C'], ax=axE, width=0.75, legend=False)
    # axE = sns.lineplot(data=df, x='SNP type', y='Percentage of SMs', hue='Layer', palette=colour, ax=axE, legend=False)
    axE.axvline(1.5, linestyle='dashed', color='grey', linewidth=0.75)
    axE.set_ylim([0, 100])
    axE.set_xlim([-0.5, 5.5])
    axE = pretty_ticks(0, 100, 5, axE, 'y')

    axE.set_ylabel('SMs (in %)')
    axE.set_xlabel('')
    axE.tick_params(axis='x', rotation=45)
    axE = cleanax(axE)
    axE.plot([], [], linewidth=5, color=colour['L1', 'SNV'], label='L1')
    axE.plot([], [], linewidth=5, color=colour['L2', 'SNV'], label='L2')
    axE.plot([], [], linewidth=5, color='darkgrey', label='SNV')
    axE.plot([], [], linewidth=5, color='lightgrey', label='Indel')
    axE.legend(frameon=False, bbox_to_anchor=[0.5, 0.2])
    axE.text(0, 1, 'D', transform=axE.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    # </editor-fold>


    # <editor-fold desc="D: Get distibution of SMs SNPs in genomic triplets. Are CpGs enriched?">
    snppos = datafilter.loc[datafilter.type == 'SNV'].drop_duplicates(subset='chromosome position alt_allele'.split())

    # DO NOT DELETE
    # # Count the total number of k-mers
    # threemers = sorted(canonical_kmers(3))
    # args_list = [Namespace(fasta=Namespace(name=refcur), kmer=kmer, canonical=True, loc=False) for kmer in threemers]
    # with Pool(processes=4) as pool:
    #     kmercnts = pool.map(cntkmer, args_list)
    # kmercnts = dict(zip(threemers, kmercnts))
    #
    # calculated kmer values
    kmercnts = {'AAA': 19620375,
                'AAC': 8618744,
                'AAG': 9552536,
                'AAT': 14242646,
                'ACA': 8596065,
                'ACC': 5517024,
                'ACG': 2579812,
                'ACT': 7163269,
                'AGA': 9019650,
                'AGC': 5113096,
                'AGG': 5649131,
                'ATA': 11234596,
                'ATC': 7388202,
                'ATG': 9340964,
                'CAA': 11953127,
                'CAC': 5617513,
                'CAG': 5065838,
                'CCA': 7975032,
                'CCC': 4216506,
                'CCG': 2013804,
                'CGA': 3194618,
                'CGC': 1528087,
                'CTA': 5900178,
                'CTC': 6426595,
                'GAA': 10101532,
                'GAC': 4260053,
                'GCA': 5929901,
                'GCC': 3635269,
                'GGA': 6485670,
                'GTA': 5359856,
                'TAA': 10359257,
                'TCA': 9476452}

    refseq = readfasta(refcur)
    trip = deque()
    for p in snppos.itertuples(index=False):
        trip.append(refseq[p.chromosome][p.position-2: p.position+1])
    snppos['trip_reg'] = trip

    # Testing for layers separately didn't return in any statistically significant 3mers. Reason: Too few SNVs
    # So, checking the counts for the combined table
    # stats = deque()
    # s = snppos.copy()
    # snpcnt = s.shape[0]
    # for k in kmercnts:
    #     cnt = s.loc[s.trip_reg == k].shape[0] + s.loc[s.trip_reg == revcomp(k)].shape[0]
    #     stats.append((l, k, *fisher_exact(np.array([[cnt, snpcnt], [kmercnts[k], 233151596]]), alternative='greater')))
    # stats = pd.DataFrame(stats)
    # stats.columns = 'Layer kmer odds pvalue'.split()
    # stats['p.val.adjusted'] = multipletests(stats.pvalue.to_numpy(), method='fdr_bh')[1]
    # print(f"Kmers enriched for somatic mutations")
    # print(stats.loc[stats['p.val.adjusted'] < 0.05])
    # This shows that ACG, AGG, CGA are significantly higher

    # Compare L1 vs L2 distribution
    stats = deque()
    for l in 'L1 L2'.split():
        s = snppos.loc[snppos.Layer == l].copy()
        snpcnt = s.shape[0]
        stats.append([(s.loc[s.trip_reg == k].shape[0] + s.loc[s.trip_reg == revcomp(k)].shape[0])/snpcnt for k in kmercnts])
    print("Triplet distribution between L1 and L2 is not statistically different.")
    print(ttest_rel(stats[0], stats[1]))

    tripcntabs = defaultdict()
    tripcntnorm = defaultdict()
    # for l in 'L1 L2 shared'.split():
    for l in 'L1 L2'.split():
        trip = snppos.loc[snppos.Layer == l, 'trip_reg']
        tripcnt = Counter(trip)
        tripkeys = list(tripcnt.keys())
        tripcntclean = defaultdict(int)
        for k in tripkeys:
            kr = revcomp(k)
            kmin = min(k, kr)
            if kmin in tripcntclean:
                continue
            tripcntclean[kmin] = 0
            try:
                tripcntclean[kmin] += tripcnt[k]
                tripcnt.pop(k)
            except KeyError:
                pass
            try:
                tripcntclean[kmin] += tripcnt[kr]
                tripcnt.pop(kr)
            except KeyError:
                pass
        tripcntabs[l] = {k: tripcntclean[k] for k in kmercnts}
        tripcntnorm[l] = {k: tripcntclean[k]/kmercnts[k] for k in kmercnts}
    cpgkeys = set([k for v in tripcntnorm.values() for k in v.keys() if 'CG' in k])
    tripletjson = defaultdict()
    tripletjson['triplet counts'] = tripcntabs
    tripletjson['triplet count normalised'] = {k: {k1: np.format_float_scientific(v1, unique=False, precision=3) for k1, v1 in v.items()} for k, v in tripcntnorm.items()}
    # tripletjson['cgsnpcnt'] = {l: sum([tripcntabs[l][c] for c in cpgkeys]) for l in 'L1 L2 shared'.split()}
    tripletjson['cgsnpcnt'] = {l: sum([tripcntabs[l][c] for c in cpgkeys]) for l in 'L1 L2'.split()}
    # tripletjson['othersnpcnt'] = {l: sum([v for k, v in tripcntabs[l].items() if k not in cpgkeys]) for l in 'L1 L2 shared'.split()}
    tripletjson['othersnpcnt'] = {l: sum([v for k, v in tripcntabs[l].items() if k not in cpgkeys]) for l in 'L1 L2'.split()}
    with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/axillary_plots/triplet_sm_stats.json', 'w') as fout:
        json.dump(tripletjson, fout, indent=2)
    # tripkeys = sorted(kmercnts, key=lambda x: sum([tripcntnorm[l][x] for l in 'L1 L2 shared'.split()]), reverse=True)
    tripkeys = sorted(kmercnts, key=lambda x: sum([tripcntnorm[l][x] for l in 'L1 L2'.split()]), reverse=True)
    tripkeys = [t for t in tripkeys if tripcntnorm['L1'][t] + tripcntnorm['L2'][t] > 0]
    axwidth = (barwidth*len(tripkeys)/pltwidth)/1
    axD = fig1.add_axes([0.5, 0.7, axwidth, 0.1])
    ybottom = np.zeros_like(tripkeys, dtype='float')
    for l in 'L1 L2'.split():
        y = [tripcntnorm[l][k] for k in tripkeys]
        axD.bar(tripkeys, y, bottom=ybottom, width=0.75, label=l, color=colour[l])
        ybottom += y
    axD.ticklabel_format(axis='y', useOffset=False, style='plain')
    axD.set_yticklabels(axD.get_yticks() * 1000000)
    axD.set_xlim([-0.5, len(tripkeys)-0.5])
    axD.set_xlabel('Triplet context')
    axD.set_ylabel('SM ratio (x$10^{-6}$)')
    axD = cleanax(axD)
    axD.tick_params(axis='x', rotation=45)
    axD.text(0, 1, 'E', transform=axD.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')

    # </editor-fold>


    # <editor-fold desc="E: Indel size distribution">
    inddf = datafilter.loc[datafilter['type'] == 'Indel'].drop_duplicates(subset=['chromosome', 'position', 'ref_allele', 'alt_allele', 'Layer'])
    indsize = dict()
    for grp in inddf.groupby('Layer'):
        indsize[grp[0]] = defaultdict(int, Counter([len(a)-1 if '+' in a else 1-len(a) for a in grp[1].alt_allele]))
    maxx = max(map(abs, unlist([list(v.keys()) for v in indsize.values()])))
    xran = list(range(-maxx, maxx+1))
    allx = list(product(['L1', 'L2'], xran))
    indsizeall = [(v[0], v[1], indsize[v[0]][v[1]]) for v in allx]
    indsizeall = pd.DataFrame(indsizeall)
    indsizeall.columns = 'tissue size count'.split()
    axwidth = barwidth*10/pltwidth
    axG = fig1.add_axes([0.075, 0.5, axwidth, 0.125])
    axG = sns.lineplot(indsizeall, x='size', y='count', hue='tissue', palette=[colour[l, 'SNV'] for l in 'L1 L2'.split()])
    axG.set_xticks([-10, -2, 0, 2 , 10])
    axG.xaxis.get_majorticklabels()[1].set_horizontalalignment('right')
    axG.xaxis.get_majorticklabels()[3].set_horizontalalignment('left')
    axG.legend().remove()
    axG = cleanax(axG)
    axG.axvline(2, linestyle='dashed', color='grey', linewidth=0.75, zorder=0)
    axG.axvline(-2, linestyle='dashed', color='grey', linewidth=0.75, zorder=0)
    axG.set_xlabel("Size (in bp)")
    axG.set_ylabel("Indels")
    axG.text(0, 1, 'F', transform=axG.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')

    # </editor-fold>


    # <editor-fold desc="F: SMs overlapping genes, TEs, and intergenic, annotations">

    # # <editor-fold desc="DO NOT DELETE. Get uncallable regions. Ran once not needed again.">
    # def get_uncallable_regions(fi, fo, rng):
    #     with open(fi, 'r') as fin, open(fo, 'w') as fout:
    #         for line in tqdm(fin):
    #             line = line.strip().split()
    #             if rng[0] <= int(line[3]) <= rng[1]: continue
    #             fout.write(f'{line[0]}\t{int(line[1])-1}\t{int(line[1])}\n')
    #     return
    # # END
    # # Get positions that are uncallable in each of the layer samples
    # fins2 = deque()
    # for n, b in bnamedict.items():
    #     for l in 'l1 l2'.split():
    #         fins2.append([f'./layer_samples/{n}/{n}_{l}/bam_read_counts_b30_q10.bt2.txt', f"{b}_{l}_uncallable_positions.bed", scovs[n][l]])
    # with Pool(processes=7) as pool:
    #     pool.starmap(get_uncallable_regions, fins2)
    #
    # # Merge the BED regions
    # for n, b in bnamedict.items():
    #     for l in 'l1 l2'.split():
    #         with open(f'{b}_{l}_uncallable_positions.merged.bed', 'w') as fout:
    #             p = Popen(f'bedtools merge -i {b}_{l}_uncallable_positions.bed '.split(), stdout=fout)
    #             p.wait()
    # # </editor-fold>

    gff = bt.BedTool('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.pasa_out.sort.protein_coding.3utr.gff3').saveas()
    te_rep = bt.BedTool('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/repeat/RepeatMasker/cur.genome.v1.fasta.ann.gff3').saveas()

    genes = gff.filter(lambda x: x[2] == 'gene').saveas()
    cds = gff.filter(lambda x: x[2] == 'CDS').sort().saveas()
    intron = genes.subtract(gff.filter(lambda x: x[2] not in {'gene', 'mRNA'})).sort().saveas()
    utr = gff.filter(lambda x: x[2] in {'five_prime_UTR', 'three_prime_UTR'}).sort().saveas()
    # Get intergenic as the genomic region other than CDS + UTR + INTRON + TE_REP
    intergenic = bt.BedTool('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.bed').filter(lambda x: 'CUR' in x[0]).saveas()
    intergenic = intergenic.subtract(cds.cat(utr).cat(intron).cat(te_rep).sort().merge())

    pos = datafilter.copy()
    pos['start'] = pos.position - 1
    pos['end'] = pos.position
    pos = pos[['chromosome', 'start', 'end', 'alt_allele', 'Layer', 'type']].copy()
    pos.drop_duplicates(inplace=True)
    posbt = bt.BedTool.from_dataframe(pos)
    genepos = posbt.intersect(genes, u=True).to_dataframe()
    cdspos = posbt.intersect(cds, u=True).to_dataframe()
    cdspos.columns = pos.columns
    cdspos['anno'] = 'cds'
    intronpos = posbt.intersect(intron, u=True).to_dataframe()
    intronpos.columns = pos.columns
    intronpos['anno'] = 'intron'
    utrpos = posbt.intersect(utr, u=True).to_dataframe()
    utrpos.columns = pos.columns
    utrpos['anno'] = 'utr'
    tepos = posbt.intersect(te_rep, u=True).to_dataframe()
    tepos.columns = pos.columns
    tepos['anno'] = 'te'

    anno = np.array(['intergenic']*pos.shape[0])
    g = pos.merge(tepos, how='left', on=['chromosome', 'start', 'end', 'alt_allele', 'Layer', 'type'])
    anno[g.anno == 'te'] = 'te'
    g = pos.merge(intronpos, how='left', on=['chromosome', 'start', 'end', 'alt_allele', 'Layer', 'type'])
    anno[g.anno == 'intron'] = 'intron'
    g = pos.merge(utrpos, how='left', on=['chromosome', 'start', 'end', 'alt_allele', 'Layer', 'type'])
    anno[g.anno == 'utr'] = 'utr'
    g = pos.merge(cdspos, how='left', on=['chromosome', 'start', 'end', 'alt_allele', 'Layer', 'type'])
    anno[g.anno == 'cds'] = 'cds'
    pos['anno'] = anno

    annotypes = 'cds utr intron te intergenic'.split()
    regsize = defaultdict(deque)
    for n, b in bnamedict.items():
        print(b)
        for l in 'l1 l2'.split():
            regsize['cds'].append(sum([int(b[2])-int(b[1])+1 for b in cds.sort().merge().subtract(f'{b}_{l}_uncallable_positions.merged.bed')]))
            regsize['utr'].append(sum([int(b[2])-int(b[1])+1 for b in utr.sort().merge().subtract(f'{b}_{l}_uncallable_positions.merged.bed')]))
            regsize['intron'].append(sum([int(b[2])-int(b[1])+1 for b in intron.sort().merge().subtract(f'{b}_{l}_uncallable_positions.merged.bed')]))
            regsize['te'].append(sum([int(b[2])-int(b[1])+1 for b in te_rep.sort().merge().subtract(f'{b}_{l}_uncallable_positions.merged.bed')]))
            regsize['intergenic'].append(sum([int(b[2])-int(b[1])+1 for b in intergenic.sort().merge().subtract(f'{b}_{l}_uncallable_positions.merged.bed')]))


    # Fisher's exact test shows that the CDS have lower than expected number of SMs
    stats = deque()
    for l in 'L1 L2'.split():
        poscnt = pos.loc[pos.Layer == l].anno.value_counts()
        for a in annotypes:
            stats.append((l, a, poscnt[a], np.mean(regsize[a]), *fisher_exact(np.array([[poscnt[a], np.mean(regsize[a])], [sum(poscnt), chrom_genome_size]]), alternative='greater')))
    stats = pd.DataFrame(stats)
    stats.columns = 'Layer anno count length odds pvalue'.split()
    stats['padjusted'] = multipletests(stats.pvalue, method='fdr_bh')[1]
    print(stats.loc[stats.padjusted < 0.05])

    xticks = ['cds', 'utr', 'intron', 'te', 'intergenic']
    axwidth = barwidth*10/pltwidth
    axH = fig1.add_axes([0.075, 0.3, axwidth, 0.125])
    bar_width = 0.75/2
    for i, l in enumerate('L1 L2'.split()):
        y_bottom = [0]*5
        for t in 'SNV Indel'.split():
            cnts = pos.loc[(pos.Layer == l) & (pos.type == t), 'anno'].value_counts()
            y = [(cnts[x]*1000000)/np.mean(regsize[x]) if x in cnts else 0 for x in xticks]
            axH.bar(np.arange(5) + bar_width*(i - 0.5), y, label=f'{l} {t}', bottom=y_bottom, color=colour[l, t], width=bar_width)
            y_bottom = y
    axH.set_xticks(np.arange(5))
    axH.set_xticklabels(['Coding', 'UTR', 'Intron', 'TE/Repeat', 'Intergenic'])
    axH.set_ylabel('SMs (per MBp)')
    axH.set_xlabel('')
    axH = cleanax(axH)
    axH.set_xlim(-0.5, 4.5)
    axH.tick_params(axis='x', rotation=45)
    axH.text(0, 1, 'G', transform=axH.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')


    fig2 = plt.figure(figsize=[4, 4], dpi=300)
    ax = fig2.add_subplot(2, 1, 1)
    for i, l in enumerate('L1 L2'.split()):
        y_bottom = [0]*5
        for t in 'SNV Indel'.split():
            cnts = pos.loc[(pos.Layer == l) & (pos.type == t), 'anno'].value_counts()
            y = [cnts[x] if x in cnts else 0 for x in xticks]
            ax.bar(np.arange(5) - (bar_width*(1-i)), y, label=f'{l} {t}', bottom=y_bottom, color=colour[l, t], width=bar_width)
            y_bottom = y
        ax.set_ylabel('Mutation count')
        ax.set_xlabel('')
        ax.set_xticks(np.arange(5))
        ax.set_xticklabels(['Coding', 'UTR', 'Intron', 'TE/Repeat', 'Intergenic'])
        ax.legend(frameon=False, loc=2)
        ax = cleanax(ax)
    ax = fig2.add_subplot(2, 1, 2)
    for i, l in enumerate('L1 L2'.split()):
        y_bottom = [0]*5
        for t in 'SNV Indel'.split():
            cnts = pos.loc[(pos.Layer == l) & (pos.type == t), 'anno'].value_counts()
            cnts = (cnts/pos.Layer.value_counts()[l])*100
            y = [cnts[x] if x in cnts else 0 for x in xticks]
            ax.bar(np.arange(5) - (bar_width*(1-i)), y, label=f'{l} {t}', bottom=y_bottom, color=colour[l, t], width=bar_width)
            y_bottom = y
            ax.set_xticks(np.arange(5))
            ax.set_xticklabels(['Coding', 'UTR', 'Intron', 'TE/Repeat', 'Intergenic'])
            ax.set_ylabel('Percent of mutations')
            ax.set_xlabel('Annotation')
            ax = cleanax(ax)
    plt.tight_layout(pad=0)
    plt.savefig(f'{cwd}/../all_layer_somatic_variants.annotation_overlap.supplementary.png')
    plt.savefig(f'{cwd}/../all_layer_somatic_variants.annotation_overlap.supplementary.pdf')
    plt.close()
    # subprocess.call(f'inkscape {cwd}/all_layer_somatic_variants.annotation_overlap.svg -M {cwd}/all_layer_somatic_variants.annotation_overlap.emf', shell=True)
    # </editor-fold>


    # <editor-fold desc="G: Get syri plot with SMs">
    df = datafilter.drop_duplicates(subset=['chromosome', 'position', 'ref_allele', 'alt_allele'])
    for grp in df.groupby('Layer'):
        pos = grp[1].copy()
        pos = pos['chromosome position'.split()]
        pos.columns = 'chromosome end'.split()
        pos['start'] = pos.end - 1
        pos = pos['chromosome start end'.split()]
        pos.to_csv(f'{cwd}/{grp[0]}_sm_pos.bed', index=False, header=False, sep='\t')
    os.chdir(cwd)
    args = Namespace(H=None, R=False, S=0.4, W=None, b='agg', bp=None, cfg=Namespace(name=f'{cwd}/../base.cfg'), chr=None, chrname=Namespace(name=f'{cwd}/../plotsr_chrnames.txt'), chrord=None, d=300, f=SMALL_SIZE, genomes=Namespace(name=f'{cwd}/../genomes.txt'), itx=False, log='WARN', logfin=Namespace(name='plotsr.log'), markers=None, nodup=False, noinv=False, nosyn=False, notr=False, o='plotsr.pdf', reg=None, rtr=False, s=25000, sr=[Namespace(name='syri.out')], tracks=Namespace(name='tracks.txt'), v=True)
    axI = fig1.add_axes([0.35, 0.3, 0.65, 0.225])
    axI = inlineplotsr(axI, args)
    axI = cleanax(axI)
    axI.spines[:].set_visible(False)
    axI.text(0, 1, 'H', transform=axI.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    axI.set_xlabel('Chromosome ID')
    ## Run plotsr using the above BED files as tracks (/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run)
    ## plotsr command:
    ## plotsr --sr syri.out --genomes genomes.txt --tracks tracks.txt -o syri_sm_tracks.png -S 0.3 -R -W 6 -H 8 -s 25000 -f 10
    # </editor-fold>


    # <editor-fold desc="H: Get syntenic and SR regions from syri.out for sm_overlap_with_structural_annotation">

    curanno = deque()
    with open(syriout, 'r') as f:
        for line in f:
            line = line.strip().split()
            if line[10] in {'SYN', 'INV', 'TRANS', 'INVTR', 'DUP', 'INVDP'}:
                curanno.append(line[:3] + [line[10]])
                continue
            if line[10] == 'NOTAL':
                if line[0] != '-':
                    curanno.append(line[:3] + [line[10]])
                # break
    curanno = pd.DataFrame(curanno)
    curanno.columns = 'chromosome start end anno'.split()
    curanno.loc[:, 'start end'.split()] = curanno.loc[:, 'start end'.split()].astype(int)
    cursyn = curanno.loc[curanno['anno'] == 'SYN']
    cursyn = bt.BedTool.from_dataframe(cursyn)
    cursr = curanno.loc[curanno['anno'] != 'SYN']
    cursr = bt.BedTool.from_dataframe(cursr)
    cursynlen = 0
    cursrlen = 0
    for grp in curanno.groupby('chromosome'):
        cursynlen += sumranges(mergeRanges(np.array(grp[1].loc[grp[1]['anno'] == 'SYN']['start end'.split()])))
        cursrlen += sumranges(mergeRanges(np.array(grp[1].loc[grp[1]['anno'] != 'SYN']['start end'.split()])))

    smdist = deque()
    # for tissue in 'L1 L2 shared'.split():
    for tissue in 'L1 L2'.split():
        sms = bt.BedTool(f'{cwd}/{tissue}_sm_pos.bed').saveas()
        sms = sms.filter(func=lambda x: 'cur' not in x[0]).saveas()
        synpos = sms.intersect(cursyn, u=True)
        srpos = sms.intersect(cursr, u=True)
        smdist.append([tissue, round(len(synpos)/len(sms), 2), round(len(srpos)/len(sms), 2), np.format_float_scientific((len(synpos)/len(sms))/cursynlen, unique=False, precision=3), np.format_float_scientific((len(srpos)/len(sms))/cursrlen, unique=False, precision=3)])
    smdist = pd.DataFrame(smdist)
    smdist = smdist[[0, 3, 4]]
    smdist.columns = 'tissue syntenic rearranged'.split()
    smdist = pd.melt(smdist, id_vars='tissue')
    smdist.value = smdist.value.astype('float')

    axwidth = barwidth*4/pltwidth
    bar_width = 0.75/2
    axJ = fig1.add_axes([0.9, 0.5, axwidth, 0.1])
    axJ = sns.barplot(data=smdist, x='variable', y='value', hue='tissue', hue_order='L1 L2'.split(), palette=colour, ax=axJ)
    axJ.ticklabel_format(axis='y', useOffset=False, style='plain')
    axJ.get_legend().remove()
    axJ.set_ylim([0, 0.7e-8])
    axJ.set_yticks(np.arange(0, 0.8, 0.25)/100000000)
    axJ.set_yticklabels(np.arange(0, 0.8, 0.25))
    axJ.set_ylabel('SMs (x$10^{-8}$/Mbp)')
    axJ.set_xlabel('')
    axJ.tick_params(axis='x', rotation=45)
    axJ = cleanax(axJ)
    axJ.text(0, 1, 'I', transform=axJ.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    # </editor-fold>

    os.chdir('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/')


    plt.savefig(f"fig1_v2.png", dpi=300)
    plt.savefig(f"fig1_v2.pdf", dpi=300)
    plt.close()
    # </editor-fold>


    # <editor-fold desc="OBSOLETE: Figure 2 plots">
    barwidth = 0.125
    pltwidth = 4.76
    dpi = 100
    fig2 = plt.figure(figsize=[4.76, 6], dpi=dpi)
    trans = mtransforms.ScaledTranslation(-5/dpi, 5/dpi, fig2.dpi_scale_trans)

    # <editor-fold desc="A: Tree plot">

    tree = mpimg.imread(f'{cwd}/../tree_sketch_2.png')
    axA = fig2.add_axes([0.03, 0.575, 0.608, 0.4])
    axA.imshow(tree, aspect='auto')
    axA.text(0, 1, 'A', transform=axA.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    axA.spines[:].set_visible(False)
    axA.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)

    # </editor-fold>


    # <editor-fold desc="B: Dendrogram for SMs in L1 and L2">
    # hierarchical clustering between branches (this would show whether the SMs fit the morphology of the tree)
    # fig = plt.figure(figsize=[3, 5], dpi=300)
    axB1 = fig2.add_axes([0.675, 0.77, 0.325, 0.18])
    axB2 = fig2.add_axes([0.675, 0.45, 0.325, 0.18])
    axdict = {0: axB1, 1: axB2}
    for i, l in enumerate('L1 L2'.split()):
        pos = datafilter.loc[datafilter.Layer.isin([l, 'shared']), ['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele']].copy()
        # pos = datafilter.loc[datafilter.Layer == l, ['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele']].copy()
        alts = {(row[0], row[1], row[3]) for row in pos.itertuples(index=False)}
        branchpos = np.zeros([7, len(alts)], dtype=int)

        for j, branch in enumerate(branches):
            bdf = pos.loc[pos.branch == branch]
            balts = {(row[0], row[1], row[3]) for row in bdf.itertuples(index=False)}
            for k, alt in enumerate(alts):
                if alt in balts:
                    branchpos[j, k] = 1
        ## Get distances for hierarchical clustering
        AM = pdist(branchpos, metric='hamming')
        # AM = pdist(branchpos, lambda u, v : 1 - (sum((u==1) & (v==1))/len(u)))
        Z = linkage(AM, method='ward')
        ax = axdict[i]
        # dendrogram(optimal_leaf_ordering(Z, AM), link_color_func=lambda x: 'black', ax=ax, leaf_font_size=SMALL_SIZE)
        garb = dendrogram(Z, link_color_func=lambda x: 'black', ax=ax, leaf_font_size=SMALL_SIZE)
        ax.spines[:].set_visible(False)
        ax.tick_params(left=False, labelleft=False)
        # ax.set_xticks(labels=[branches[i] for i in leaves_list(optimal_leaf_ordering(Z, AM))], ticks=ax.get_xticks(), rotation=90)
        ax.set_xticks(labels=[branches[i] for i in garb['leaves']], ticks=ax.get_xticks(), rotation=90)
        ax.set_title(l)
    axB1.text(0, 1, 'B', transform=axB1.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')

    # # Get maximum spanning tree
    # # Get adjacency matrix for minimum spanning tree
    # AM2 = 30 - squareform(AM)
    # np.fill_diagonal(AM2, 0)
    # g = ig.Graph.Weighted_Adjacency(AM2, mode="undirected")
    # g.vs['name'] = branches
    # inv_weight = [1./w for w in g.es["weight"]]
    # T = g.spanning_tree(weights=inv_weight)
    # print(T.is_connected())
    #
    # ax = fig.add_subplot(2, 2, 2)
    # ig.plot(g, target=ax, layout='davidson_harel', vertex_label=g.vs["name"], edge_width=np.log1p(g.es['weight']), vertex_size=0.2)
    # ax.set_xlabel('branch connectivity')
    # ax = fig.add_subplot(2, 2, 4)
    # ig.plot(T, target=ax, layout='davidson_harel', vertex_label=T.vs["name"], edge_width=np.log1p(T.es['weight']), vertex_size=0.2)
    # ax.set_xlabel('Maximum spanning tree')
    # plt.tight_layout(h_pad=2, w_pad=3)
    # plt.savefig(f'{cwd}/sm_branch_clustering.png')
    # plt.close()

    # </editor-fold>


    # <editor-fold desc="C: Grouping of SMs based on branching">
    # Generate the Figure 5 from Tomimoto and Satake, 2023
    l1cnts = defaultdict(int)
    l2cnts = defaultdict(int)
    shcnts = defaultdict(int)
    for grp in datafilter.groupby('chromosome position alt_allele Layer'.split()):
        selb = tuple([1 if branch in grp[1].branch.to_list() else 0 for branch in branchorder])
        if grp[0][3] == 'L1':
            l1cnts[selb] += 1
        elif grp[0][3] == 'L2':
            l2cnts[selb] += 1
        elif grp[0][3] == 'shared':
            shcnts[selb] += 1

    bkeys = set(list(l1cnts) + list(l2cnts))
    bkeys = sorted(sorted(bkeys), key=lambda x: sum(x))

    cnts = {'L1': l1cnts, 'L2': l2cnts, 'shared': shcnts}
    sharedcnts = dict()
    sharedcnts['L1'] = sum([v for k, v in cnts['L1'].items() if 1 < sum(k) < 7])
    sharedcnts['L2'] = sum([v for k, v in cnts['L2'].items() if 1 < sum(k) < 7])

    axwidth = barwidth*14/pltwidth
    axC1 = fig2.add_axes([0.12, 0.35, axwidth, 0.125])
    axC2 = fig2.add_axes([0.12, 0.2, axwidth, 0.125])
    axC3 = fig2.add_axes([0.12, 0.05, axwidth, 0.125])
    axdict = {0: axC1, 1: axC2}
    for i, l in enumerate('L1 L2'.split()):
        data = cnts[l]
        total_cnt = sum(data.values())
        ax = axdict[i]
        ax.bar(x=[str(b) for b in bkeys], height=[data[k]/total_cnt for k in bkeys], color=colour[l], width=0.75)
        ax.tick_params(axis='x',
                       which='both',
                       bottom=False,
                       top=False,
                       labelbottom=False)
        ax.set_ylim([0, 0.25])
        ax.set_ylabel(f'{l} SMs (in %)')
        ax.vlines(x=6.5, ymin=0, ymax=25, color='black', linestyle='dotted')
        ax=cleanax(ax)
        ax.grid(which='major', axis='x')
        ax.set_axisbelow(True)
        ax.set_xlim([-0.5, 13.5])

    xpos, ypos = np.where(np.array(bkeys)==1)
    ypos = 6 - ypos
    axC3.scatter([x for i, x in enumerate(xpos) if x not in {8, 12}], [y for i, y in enumerate(ypos) if xpos[i] not in {8, 12}], marker=r'$\checkmark$', color='black')
    axC3.scatter([x for i, x in enumerate(xpos) if x in {8, 12}], [y for i, y in enumerate(ypos) if xpos[i] in {8, 12}], marker=r'$\chi$', color='red')
    axC3.set_yticks(ticks=range(7))
    axC3.set_xticks(ticks=range(14))
    axC3.set_yticklabels(labels=branchorder[::-1])
    axC3.grid()
    axC3.set_axisbelow(True)
    axC3.set_xlim([-0.5, 13.5])
    axC3.set_ylim([-0.5, 6.5])
    axC3.tick_params(axis='x',
                     which='both',
                     bottom=False,
                     top=False,
                     labelbottom=False)
    axC3.set_xlabel('Selected branches')
    axC3.spines[:].set_visible(False)
    axC1.text(0, 1, 'C', transform=axC1.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')

    # </editor-fold>


    # <editor-fold desc="D: Correlation of SM count against branch length and branching events count from the primary branching">
    branchlength = {'wt_1': (0.93 + 1.25 + 1.5),
                    'wt_7': (0.93 + 1.25 + 1.0),
                    'wt_18': (0.93 + 2.25),
                    'wt_19': (3.13),
                    'mut_11_1': (0.77 + 0.92 + 1.51 + 0.1),
                    'mut_11_2': (0.77 + 0.92 + 1.51 + 0.2),
                    'mut_15': (0.77 + 0.92 + 0.08 + 1.45)}
    # branchcount = {'wt_1': 3,
    #                'wt_7': 3,
    #                'wt_18': 3,
    #                'wt_19': 2,
    #                'mut_11_1': 4,
    #                'mut_11_2': 4,
    #                'mut_15': 4}
    df = datafilter.copy()
    df.set_index('chromosome position alt_allele'.split(), inplace=True)
    # Filter positions that are in branches

    grp = datafilter.groupby('chromosome position alt_allele'.split()).nunique()
    grp = grp.loc[grp.branch == 7]
    df = df.loc[~df.index.isin(grp.index)]
    df = df.reset_index()

    # branchsmcnt = dict()
    layersmcount = defaultdict(dict)
    for l in 'L1 L2'.split():
        layersmcount[l] = {branch: df.loc[df.branch==branch].loc[df.Layer.isin([l, 'shared'])].shape[0] for branch in branches}
    #
    # mutsinallL1 = allsmpivot.loc[:, [b+'_L1' for b in branches]]              # Get all positions that are in L1
    # allsmmatfilt = allsmpivot.loc[mutsinallL1.apply(sum, axis=1) != 7]        # Filter positions that are in all L1
    #
    # for branch in branches:
    #     bdf = allsmmatfilt.loc[:, [f'{branch}_{l}' for l in ['L1', 'L2', 'leaf']]]
    #     bdf = bdf.loc[(bdf != 0).any(axis=1)]
    #     branchsmcnt[branch] = bdf.shape[0]
    #     layersmcount['L1'][branch] = sum(allsmmatfilt[f'{branch}_L1'])
    #     layersmcount['L2'][branch] = sum(allsmmatfilt[f'{branch}_L2'])
    #     layersmcount['leaf'][branch] = sum(allsmmatfilt[f'{branch}_leaf'])

    # fig = plt.figure(figsize=[3, 3], dpi=300)
    # ax = fig.add_subplot(4, 2, 1)
    # ax = sns.regplot(x=[branchlength[b] for b in branches], y=[branchsmcnt[b] for b in branches], ax=ax, label='All', color=colour['any'])
    # r, p = sp.stats.spearmanr([branchlength[b] for b in branches], [branchsmcnt[b] for b in branches])
    # ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p), transform=ax.transAxes)
    # # ax.set_xlabel('branch length (in m)')
    # ax.set_ylabel('Total number of SM')
    # ax = cleanax(ax)
    # ax = fig.add_subplot(4, 2, 2)
    # ax = sns.regplot(x=[branchcount[b] for b in branches], y=[branchsmcnt[b] for b in branches], ax=ax, label='All', color=colour['any'])
    # r, p = sp.stats.spearmanr([branchcount[b] for b in branches], [branchsmcnt[b] for b in branches])
    # ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p), transform=ax.transAxes)
    # # ax.set_xlabel('Number of branching event')
    # ax = cleanax(ax)
    # for l in ['L1', 'L2', 'leaf']:
    # ax = fig.add_subplot(2, 1, i)
    # ax = fig.add_subplot()
    axD = fig2.add_axes([0.675, 0.075, 0.3, 0.2])
    for i, l in enumerate(['L1', 'L2']):
        axD = sns.regplot(x=[branchlength[b] for b in branches], y=[layersmcount[l][b] for b in branches], ax=axD, label=l, color=colour[l], scatter_kws={'s': 5}, ci=None)
        r, p = spearmanr([branchlength[b] for b in branches], [layersmcount[l][b] for b in branches])
        axD.text(.6, .2 - (i/10), 'r={:.2f}, p={:.2g}'.format(r, p), transform=axD.transAxes, label=l, color=colour[l])
        axD.set_ylabel(f'SMs')
        axD = cleanax(axD)
        i += 1
        # ax = fig.add_subplot(4, 2, i)
        # ax = sns.regplot(x=[branchcount[b] for b in branches], y=[layersmcount[l][b] for b in branches], ax=ax, label=l, color=colour[l])
        # r, p = sp.stats.spearmanr([branchcount[b] for b in branches], [layersmcount[l][b] for b in branches])
        # ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p), transform=ax.transAxes)
        # ax = cleanax(ax)
        # i+=1
        # # ax.set_ylabel(f'Number of SM in {l}')
        # if l == 'leaf':
        #     ax.set_xlabel('Number of branching event')
    # plt.tight_layout(h_pad=2, w_pad=2)
    axD.legend(frameon=True)
    axD.set_xlabel('Branch length (in m)')
    axD.text(0, 1, 'D', transform=axD.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    # </editor-fold>

    plt.savefig("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/fig2.png", dpi=300)
    plt.close()
    # </editor-fold>


    # <editor-fold desc="Figure 2 plots">
    barwidth = 0.125
    dpi = 100
    pltwidth = 7.24
    fig2 = plt.figure(figsize=[pltwidth, 8.5], dpi=dpi)
    trans = mtransforms.ScaledTranslation(-20/dpi, 5/dpi, fig2.dpi_scale_trans)
    constructor = DistanceTreeConstructor()
    treealgo = constructor.nj
    suptree = plt.figure()
    axC = suptree.add_subplot()

    # <editor-fold desc="A: Add tree image">
    def getcntplt(ax, branch, bar_width):
        # ax = fig1.add_axes(axispos)
        y_bottom = [0]
        ymax = 75
        ax.set_xlim([-0.5, 0.5])
        ax.set_ylim([0, ymax])
        ax.set_yticks([0, ymax])
        ax.tick_params(bottom=False, labelbottom=False)
        ax.set_xlabel(bnamedict[branch], fontdict={'backgroundcolor': 'white', 'linespacing': 0, 'horizontalalignment': 'center', 'bbox': {'pad': 1, 'fc': 'white', 'ec': 'white'}})
        ax.set_yticklabels(labels=ax.get_yticklabels(), fontdict={'backgroundcolor': 'white', 'linespacing': 0, 'horizontalalignment': 'center', 'verticalalignment': 'center', 'bbox': {'pad': 1, 'fc': 'white', 'ec': 'white'}})
        ax.spines[:].set_visible(False)
        for l in ('L1', 'L2'):
            for t in ('SNV', 'Indel'):
                y = branchscnt[(l, t)][branch]
                ax.bar([0], [y], label=f'{l} {t}', bottom=y_bottom, color=colour[(l, t)], width=bar_width)
                y_bottom = [y_bottom[0] + y]
        return ax
    # END

    tree = mpimg.imread(f'{cwd}/../tree_sketch_fruit_nolength.png')
    axA = fig2.add_axes([0.03, 0.675, 0.4, 0.3])
    axA.imshow(tree, aspect='auto')
    axA.text(0, 1, 'A', transform=axA.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    axA.spines[:].set_visible(False)
    axA.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)

    # branchscnt = {g[0]: defaultdict(int, g[1].branch.value_counts()) for g in datafilter.groupby(['Layer', 'type'])}
    # bar_width = 0.75
    # axwidth = (barwidth*1/pltwidth)/0.4
    # axA1 = getcntplt(axA.inset_axes([0.9, 0.7, axwidth, 0.1]), 'wt_1', bar_width)
    # axA2 = getcntplt(axA.inset_axes([0.9, 0.25, axwidth, 0.1]), 'wt_7', bar_width)
    # axA3 = getcntplt(axA.inset_axes([0.75, 0.85, axwidth, 0.1]), 'wt_18', bar_width)
    # axA4 = getcntplt(axA.inset_axes([0.45, 0.85, axwidth, 0.1]), 'wt_19', bar_width)
    # axA5 = getcntplt(axA.inset_axes([0.25, 0.85, axwidth, 0.1]), 'mut_15', bar_width)
    # axA6 = getcntplt(axA.inset_axes([0.1, 0.85, axwidth, 0.1]), 'mut_11_1', bar_width)
    # axA7 = getcntplt(axA.inset_axes([0.1, 0.45, axwidth, 0.1]), 'mut_11_2', bar_width)

    # </editor-fold>


    # <editor-fold desc="B: Cluster samples based on all SMs">
    axB = fig2.add_axes([0.5, 0.675, 0.2, 0.3])
    tissues = ['L1', 'L2']
    # btiss = ['_'.join([b, a]) for a in tissues for b in branches]
    btiss_lab = ['_'.join([bnamedict[b], a]) for a in tissues for b in branches]

    # samples = sorted(['_'.join(i) for i in product(branches, 'L1 L2'.split())], key=lambda x: x.rsplit('_', 1)[1])

    unism = datafilter.drop_duplicates(subset='chromosome position alt_allele'.split())
    smcnt = unism.shape[0]
    alts = set(unism.set_index(['chromosome', 'position', 'alt_allele']).index.values)
    branchpos = np.zeros([14, smcnt], dtype=int)
    for i, l in enumerate('L1 L2'.split()):
        for j, branch in enumerate(branches):
            print(branch, l, (i*7) + j)
            bdf = datafilter.loc[(datafilter.branch == branch) & (datafilter.Layer.isin([l, 'shared']))]
            balts = set(bdf.set_index(['chromosome', 'position', 'alt_allele']).index.values)
            for k, alt in enumerate(alts):
                if alt in balts:
                    branchpos[(i*7) + j, k] = 1

    AM = pdist(branchpos, metric='hamming')
    # Lower triangle distance matrix (list of lists)
    mat = [r[:i + 1] for i, r in enumerate(np.tril(squareform(AM)).tolist())]

    dm = DistanceMatrix(names=btiss_lab, matrix=mat)

    tree = treealgo(dm)
    tree.ladderize()
    tree.root.clades = tree.root.clades[::-1]
    tree.root.clades[0].clades = tree.root.clades[0].clades[::-1]
    cladel1 = tree.common_ancestor({"name": "B1_L1"}, {"name": "B7_L1"})
    cladel2 = tree.common_ancestor({"name": "B4_L2"}, {"name": "B7_L2"})
    cladel1.color = colour['L1']
    cladel2.color = colour['L2']
    pdraw(tree,
          label_func=lambda x: x if 'Inner' not in x.name else '',
          label_colors=lambda x: (colour['L1'] if 'L1' in x.name else colour['L2']) if type(x) != str else 'black',
          axes=axB)

    # Z = linkage(AM, method='ward')
    # set_link_color_palette([colour['L2'], colour['L1'], 'black'])
    # set_link_color_palette(['black'])
    # garb = dendrogram(Z, ax=axB, leaf_font_size=SMALL_SIZE, orientation='left', above_threshold_color='black', count_sort='ascending')
    axB.spines[:].set_visible(False)
    axB.tick_params(left=False, labelleft=False)
    axB.tick_params(bottom=False, labelbottom=False)
    axB.set_xlabel('')
    axB.set_ylabel('All somatic mutations')
    # ax.set_xticks(labels=[bnamedict[branches[i]] for i in leaves_list(Z)], ticks=ax.get_xticks(), rotation=0)
    # ax.set_xticks()
    # axB.set_yticks(labels=[btiss_lab[i] for i in garb['leaves']], ticks=axB.get_yticks(), rotation=0)
    axB.text(0, 1, 'B', transform=axB.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')


    # </editor-fold>

    l1allsample = list()
    l2allsample = list()
    # # <editor-fold desc="C: Cluster samples based on all SMs except those that are present in all L1">

    # axC = fig2.add_axes([0.64, 0.675, 0.1, 0.3])
    tissues = ['L1', 'L2']
    # btiss = ['_'.join([b, a]) for a in tissues for b in branches]
    btiss_lab = ['_'.join([bnamedict[b], a]) for a in tissues for b in branches]

    # samples = sorted(['_'.join(i) for i in product(branches, 'L1 L2'.split())], key=lambda x: x.rsplit('_', 1)[1])

    unism = datafilter.drop_duplicates(subset='chromosome position alt_allele'.split())
    smcnt = unism.shape[0]
    alts = set(unism.set_index(['chromosome', 'position', 'alt_allele']).index.values)
    branchpos = np.zeros([14, smcnt], dtype=int)
    for i, l in enumerate('L1 L2'.split()):
        for j, branch in enumerate(branches):
            print(branch, l, (i * 7) + j)
            bdf = datafilter.loc[(datafilter.branch == branch) & (datafilter.Layer.isin([l, 'shared']))]
            balts = set(bdf.set_index(['chromosome', 'position', 'alt_allele']).index.values)
            for k, alt in enumerate(alts):
                if alt in balts:
                    branchpos[(i * 7) + j, k] = 1
    # Remove positions present in all L1
    # for b in branches:
    #     btissm[f'{b}_L1'] = btissm[f'{b}_L1'].difference(allL1)
    print(branchpos.shape)
    l1allsample = list(np.where(sum(branchpos[:7, :]) == 7)[0])
    l2allsample = list(np.where(sum(branchpos[7:, :]) == 7)[0])
    branchpos = np.delete(branchpos, l1allsample + l2allsample, axis=1)
    print(branchpos.shape)
    # branchpos = np.delete(branchpos, np.where(sum(branchpos[7:, :]) == 7)[0], axis=1)
    print(branchpos.shape)

    AM = pdist(branchpos, metric='hamming')
    # Lower triangle distance matrix (list of lists)
    mat = [r[:i + 1] for i, r in enumerate(np.tril(squareform(AM)).tolist())]

    dm = DistanceMatrix(names=btiss_lab, matrix=mat)
    constructor = DistanceTreeConstructor()
    tree = treealgo(dm)
    tree.ladderize()
    # tree.root.clades = tree.root.clades[::-1]
    tree.root.clades[0].clades = tree.root.clades[0].clades[::-1]
    cladel1 = tree.common_ancestor({"name": "B1_L1"}, {"name": "B7_L1"})
    cladel2 = tree.common_ancestor({"name": "B4_L2"}, {"name": "B7_L2"})
    cladel1.color = colour['L1']
    cladel2.color = colour['L2']
    pdraw(tree,
          label_func=lambda x: x if 'Inner' not in x.name else '',
          label_colors=lambda x: (colour['L1'] if 'L1' in x.name else colour['L2']) if type(x) != str else 'black',
          axes=axC)

    # Z = linkage(AM, method='ward')
    # print(Z)
    # set_link_color_palette([colour['L2'], colour['L1'], 'black'])
    # set_link_color_palette(['black'])
    #
    # garb = dendrogram(Z, ax=axC, leaf_font_size=SMALL_SIZE, orientation='left',
    #            above_threshold_color='black', count_sort='descending')
    axC.spines[:].set_visible(False)
    axC.tick_params(left=False, labelleft=False)
    axC.tick_params(bottom=False, labelbottom=False)
    axC.set_xlabel('')
    axC.set_ylabel('')
    suptree.savefig('l1_l2.sm_tree.no_all_shared.pdf')
    suptree.savefig('l1_l2.sm_tree.no_all_shared.png')
    plt.close(suptree)
    # ax.set_xticks(labels=[bnamedict[branches[i]] for i in leaves_list(Z)], ticks=ax.get_xticks(), rotation=0)
    # ax.set_xticks()
    # print(leaves_list(Z))
    # axC.set_yticks(labels=[btiss_lab[i] for i in garb['leaves']], ticks=axC.get_yticks(), rotation=0)

    # </editor-fold>


    # <editor-fold desc="D: Cluster subset of samples based on all SMs except those that are present in all L1">
    # WT samples
    axD = fig2.add_axes([0.775, 0.85, 0.2, 0.125])
    tissues = ['L1', 'L2']
    # btiss = ['_'.join([b, a]) for a in tissues for b in branches]
    wtbranches = ['wt_1', 'wt_7', 'wt_18']
    btiss_lab = ['_'.join([bnamedict[b], a]) for a in tissues for b in wtbranches]

    # samples = sorted(['_'.join(i) for i in product(branches, 'L1 L2'.split())], key=lambda x: x.rsplit('_', 1)[1])

    unism = datafilter.drop_duplicates(subset='chromosome position alt_allele'.split())
    smcnt = unism.shape[0]
    alts = set(unism.set_index(['chromosome', 'position', 'alt_allele']).index.values)
    branchpos = np.zeros([len(wtbranches)*2, smcnt], dtype=int)
    for i, l in enumerate('L1 L2'.split()):
        for j, branch in enumerate(wtbranches):
            print(branch, l, (i * len(wtbranches)) + j)
            bdf = datafilter.loc[(datafilter.branch == branch) & (datafilter.Layer.isin([l, 'shared']))]
            balts = set(bdf.set_index(['chromosome', 'position', 'alt_allele']).index.values)
            for k, alt in enumerate(alts):
                if alt in balts:
                    branchpos[(i * len(wtbranches)) + j, k] = 1
    # Remove positions present in all L1
    # for b in branches:
    #     btissm[f'{b}_L1'] = btissm[f'{b}_L1'].difference(allL1)
    print(branchpos.shape)
    branchpos = np.delete(branchpos, l1allsample + l2allsample, axis=1)
    #
    # branchpos = np.delete(branchpos, np.where(sum(branchpos[:len(wtbranches), :]) == len(wtbranches))[0], axis=1)
    # print(branchpos.shape)
    # branchpos = np.delete(branchpos, np.where(sum(branchpos[len(wtbranches):, :]) == len(wtbranches))[0], axis=1)
    print(branchpos.shape)

    AM = pdist(branchpos, metric='hamming')
    # Lower triangle distance matrix (list of lists)
    mat = [r[:i + 1] for i, r in enumerate(np.tril(squareform(AM)).tolist())]

    dm = DistanceMatrix(names=btiss_lab, matrix=mat)
    constructor = DistanceTreeConstructor()
    tree = treealgo(dm)
    tree.ladderize()
    # tree.root.clades = tree.root.clades[::-1]
    # tree.root.clades[0].clades = tree.root.clades[0].clades[::-1]
    cladel1 = tree.common_ancestor({"name": "B1_L1"}, {"name": "B2_L1"})
    cladel2 = tree.common_ancestor({"name": "B1_L2"}, {"name": "B2_L2"})
    cladel1.color = colour['L1']
    cladel2.color = colour['L2']
    pdraw(tree,
          label_func=lambda x: x if 'Inner' not in x.name else '',
          label_colors=lambda x: (colour['L1'] if 'L1' in x.name else colour['L2']) if type(x) != str else 'black',
          axes=axD)
    #
    # Z = linkage(AM, method='ward')
    # print(Z)
    # set_link_color_palette([colour['L2'], colour['L1'], 'black'])
    # set_link_color_palette(['black'])
    #
    # garb = dendrogram(Z, ax=axD, leaf_font_size=SMALL_SIZE, orientation='left',
    #            above_threshold_color='black', count_sort='descending')
    axD.spines[:].set_visible(False)
    axD.tick_params(left=False, labelleft=False)
    axD.tick_params(bottom=False, labelbottom=False)
    axD.set_xlabel('')
    axD.set_ylabel('Mutations in right sub-tree')
    # ax.set_xticks(labels=[bnamedict[branches[i]] for i in leaves_list(Z)], ticks=ax.get_xticks(), rotation=0)
    # ax.set_xticks()
    # print(leaves_list(Z))
    # axD.set_yticks(labels=[btiss_lab[i] for i in garb['leaves']], ticks=axD.get_yticks(), rotation=0)
    axD.text(0, 1, 'C', transform=axD.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')

    # </editor-fold>


    # <editor-fold desc="E: Cluster subset of samples based on all SMs except those that are present in all L1">
    # MUT samples
    axE = fig2.add_axes([0.775, 0.675, 0.2, 0.125])
    tissues = ['L1', 'L2']
    # btiss = ['_'.join([b, a]) for a in tissues for b in branches]
    mutbranches = ['mut_11_1', 'mut_11_2', 'mut_15']
    btiss_lab = ['_'.join([bnamedict[b], a]) for a in tissues for b in mutbranches]

    # samples = sorted(['_'.join(i) for i in product(branches, 'L1 L2'.split())], key=lambda x: x.rsplit('_', 1)[1])

    unism = datafilter.drop_duplicates(subset='chromosome position alt_allele'.split())
    smcnt = unism.shape[0]
    alts = set(unism.set_index(['chromosome', 'position', 'alt_allele']).index.values)
    branchpos = np.zeros([len(mutbranches)*2, smcnt], dtype=int)
    for i, l in enumerate('L1 L2'.split()):
        for j, branch in enumerate(mutbranches):
            print(branch, l, (i * len(mutbranches)) + j)
            bdf = datafilter.loc[(datafilter.branch == branch) & (datafilter.Layer.isin([l, 'shared']))]
            balts = set(bdf.set_index(['chromosome', 'position', 'alt_allele']).index.values)
            for k, alt in enumerate(alts):
                if alt in balts:
                    branchpos[(i * len(mutbranches)) + j, k] = 1
    # Remove positions present in all L1
    # for b in branches:
    #     btissm[f'{b}_L1'] = btissm[f'{b}_L1'].difference(allL1)
    print(branchpos.shape)
    branchpos = np.delete(branchpos, l1allsample + l2allsample, axis=1)
    # branchpos = np.delete(branchpos, np.where(sum(branchpos[:len(mutbranches), :]) == len(mutbranches))[0], axis=1)
    # print(branchpos.shape)
    # branchpos = np.delete(branchpos, np.where(sum(branchpos[len(mutbranches):, :]) == len(mutbranches))[0], axis=1)
    print(branchpos.shape)

    AM = pdist(branchpos, metric='hamming')
    # Lower triangle distance matrix (list of lists)
    mat = [r[:i + 1] for i, r in enumerate(np.tril(squareform(AM)).tolist())]

    dm = DistanceMatrix(names=btiss_lab, matrix=mat)
    constructor = DistanceTreeConstructor()
    tree = treealgo(dm)
    tree.ladderize()
    tree.root.clades = tree.root.clades[::-1]
    # tree.root.clades[0].clades = tree.root.clades[0].clades[::-1]
    cladel1 = tree.common_ancestor({"name": "B5_L1"}, {"name": "B7_L1"})
    cladel2 = tree.common_ancestor({"name": "B5_L2"}, {"name": "B7_L2"})
    cladel1.color = colour['L1']
    cladel2.color = colour['L2']
    pdraw(tree,
          label_func=lambda x: x if 'Inner' not in x.name else '',
          label_colors=lambda x: (colour['L1'] if 'L1' in x.name else colour['L2']) if type(x) != str else 'black',
          axes=axE)
    # Z = linkage(AM, method='ward')
    # set_link_color_palette([colour['L2'], colour['L1'], 'black'])
    # set_link_color_palette(['black'])
    #
    # garb = dendrogram(Z, ax=axE, leaf_font_size=SMALL_SIZE, orientation='left',
    #            above_threshold_color='black')
    axE.spines[:].set_visible(False)
    axE.tick_params(left=False, labelleft=False)
    axE.tick_params(bottom=False, labelbottom=False)
    axE.set_xlabel('')
    axE.set_ylabel('Mutations in left sub-tree')
    # ax.set_xticks(labels=[bnamedict[branches[i]] for i in leaves_list(Z)], ticks=ax.get_xticks(), rotation=0)
    # ax.set_xticks()
    # print(leaves_list(Z))
    # axE.set_yticks(labels=[btiss_lab[i] for i in garb['leaves']], ticks=axE.get_yticks(), rotation=0)
    axE.text(0, 1, 'D', transform=axE.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')

    # </editor-fold>


    # plt.savefig('tmp2.png')

    #
    # allsmpivot = dict()
    # for grp in allsmdata.groupby(by='chromosome position alt_allele'.split()):
    #     df = grp[1].copy()
    #     dfsamples = set(df.branch.astype(str) + '_' + df.tissue.astype(str))
    #     allsmpivot[grp[0]] = {s: (1 if s in dfsamples else 0) for s in samples}
    # allsmpivot = pd.DataFrame(allsmpivot).T
    # # isSNP = ['red' if v[2][0] in '+-' else 'blue' for v in allsmpivot.index]
    # # cmap = sns.blend_palette(["white", colour['theme1']], as_cmap=True)
    # # sns.clustermap(allsmpivot, cmap=cmap, linecolor='lightgrey', dendrogram_ratio=(0.05, 0.05), linewidths=0.01, cbar_pos=None, figsize=(8, 10), row_colors=isSNP, yticklabels=False)
    # # plt.ylabel('')
    # # plt.tight_layout()
    # # plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_heatmap_in_all_samples.png', dpi=300)
    # # plt.close()
    # #
    # # sns.clustermap(allsmpivot, cmap=cmap, linecolor='lightgrey', col_cluster=False, dendrogram_ratio=(0.05, 0.0), linewidths=0.01, cbar_pos=None, figsize=(8, 10), row_colors=isSNP, yticklabels=False)
    # # plt.ylabel('')
    # # plt.tight_layout()
    # # plt.savefig(f'{cwd}/all_sm_heatmap_in_all_samples_no_col_cluster.png', dpi=300)
    # # plt.close()
    #
    #
    # tissues = ['L1', 'L2']
    # btiss = ['_'.join([b, a]) for a in tissues for b in branches]
    # btissm = dict()
    # for c in allsmpivot.columns:
    #     btissm[c] = set((allsmpivot.loc[allsmpivot[c] == 1]).index.values)
    # allL1 = btissm[f'{branches[0]}_L1']
    # for b in branches[1:]:
    #     allL1 = allL1.intersection(btissm[f'{b}_L1'])
    # # # Remove positions present in all L1
    # # for b in branches:
    # #     btissm[f'{b}_L1'] = btissm[f'{b}_L1'].difference(allL1)
    #
    # btiss_lab = ['_'.join([bnamedict[b], a]) for a in tissues for b in branches]
    # btiss_mat = np.zeros((len(btiss), len(btiss)), dtype='int')
    # for i, a in enumerate(btiss):
    #     for j, b in enumerate(btiss):
    #         if i == j: continue
    #         btiss_mat[i][j] = len(btissm[a].intersection(btissm[b]))
    #
    # mask = np.zeros_like(btiss_mat, dtype=bool)
    # mask[np.triu_indices_from(mask)] = True
    #
    # plt.close()
    # fig = plt.figure()
    # ax = fig.add_subplot()
    # ax = sns.heatmap(btiss_mat[1:, :-1], mask=mask[1:, :-1], cmap=cmap, xticklabels=btiss_lab[:-1], yticklabels=btiss_lab[1:], linewidths=0.01, cbar=False, ax=ax)
    #
    # ax = fig.add_subplot(2, 2, 2)
    # g = ig.Graph.Weighted_Adjacency(btiss_mat, mode='undirected')
    # g.vs['bname'] = [bnamedict[i.rsplit("_", 1)[0]] for i in btiss]
    # g.vs['tname'] = [colour[i.rsplit("_", 1)[1]] for i in btiss]
    # ig.plot(g, target=ax, vertex_label=g.vs["bname"], vertex_color=g.vs["tname"], edge_width=np.log10(np.array(g.es['weight']) + 1), layout=g.layout_auto(), margin=20, bbox=[100, 100], label_size=1, vertex_size=0.5)
    # plt.tight_layout()
    # plt.savefig('tmp.png')
    #
    # sns.heatmap(btiss_mat[1:, :-1], mask=mask, cmap=cmap, xticklabels=btiss, yticklabels=btiss, dendrogram_ratio=(0.05, 0.05), cbar_pos=None, linewidths=0.01, cbar=False)
    # plt.tight_layout()
    # plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_distribution.png', dpi=300)
    # plt.close()
    #
    #
    # # axdict = {0: axB1, 1: axB2}
    # # for i, l in enumerate('L1 L2'.split()):
    # #     pos = datafilter.loc[datafilter.Layer.isin([l, 'shared']), ['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele']].copy()
    # #     # pos = datafilter.loc[datafilter.Layer == l, ['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele']].copy()
    # #     alts = {(row[0], row[1], row[3]) for row in pos.itertuples(index=False)}
    # #     branchpos = np.zeros([7, len(alts)], dtype=int)
    # #
    # #     for j, branch in enumerate(branches):
    # #         bdf = pos.loc[pos.branch == branch]
    # #         balts = {(row[0], row[1], row[3]) for row in bdf.itertuples(index=False)}
    # #         for k, alt in enumerate(alts):
    # #             if alt in balts:
    # #                 branchpos[j, k] = 1
    # #     ## Get distances for hierarchical clustering
    # #     AM = pdist(branchpos, metric='hamming')
    # #     # AM = pdist(branchpos, lambda u, v : 1 - (sum((u==1) & (v==1))/len(u)))
    # #     Z = linkage(AM, method='ward')
    # #     ax = axdict[i]
    # #     # dendrogram(optimal_leaf_ordering(Z, AM), link_color_func=lambda x: 'black', ax=ax, leaf_font_size=SMALL_SIZE)
    # #     dendrogram(Z, link_color_func=lambda x: 'black', ax=ax, leaf_font_size=SMALL_SIZE)
    # #     ax.spines[:].set_visible(False)
    # #     ax.tick_params(left=False, labelleft=False)
    # #     # ax.set_xticks(labels=[branches[i] for i in leaves_list(optimal_leaf_ordering(Z, AM))], ticks=ax.get_xticks(), rotation=90)
    # #     ax.set_xticks(labels=[bnamedict[branches[i]] for i in leaves_list(Z)], ticks=ax.get_xticks(), rotation=0)
    # #     ax.set_ylabel(l)
    #
    #
    # #
    # #
    # # g = ig.Graph.Weighted_Adjacency(btiss_mat, mode='undirected')
    # # g.vs['bname'] = [bnamedict[i.rsplit("_", 1)[0]] for i in btiss]
    # # g.vs['tname'] = [colour[i.rsplit("_", 1)[1]] for i in btiss]
    # # fig = plt.figure(figsize=[8, 8])
    # # ax = fig.add_subplot(1, 1, 1)
    # # ig.plot(g, target=ax, vertex_label=g.vs["bname"], vertex_color=g.vs["tname"], edge_width=np.log1p(g.es['weight']), layout=g.layout_auto(), margin=0, bbox=[100, 100], vertex_size=0.4)
    # # plt.tight_layout()
    # # plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_sharing.png')
    # # plt.close()
    #


    # <editor-fold desc="J: Grouping of SMs based on branching">
    # Generate the Figure 5 from Tomimoto and Satake, 2023
    l1cnts = defaultdict(int)
    l2cnts = defaultdict(int)
    shcnts = defaultdict(int)
    for grp in datafilter.groupby('chromosome position alt_allele Layer'.split()):
        selb = tuple([1 if branch in grp[1].branch.to_list() else 0 for branch in branchorder])
        if grp[0][3] == 'L1':
            l1cnts[selb] += 1
        elif grp[0][3] == 'L2':
            l2cnts[selb] += 1
        elif grp[0][3] == 'shared':
            shcnts[selb] += 1

    bkeys = set(list(l1cnts) + list(l2cnts) + list(shcnts))
    bkeys = sorted(sorted(bkeys), key=lambda x: sum(x), reverse=True)

    cnts = {'L1': l1cnts, 'L2': l2cnts, 'shared': shcnts}
    sharedcnts = dict()
    sharedcnts['L1'] = sum([v for k, v in cnts['L1'].items() if 1 < sum(k) < 7])
    sharedcnts['L2'] = sum([v for k, v in cnts['L2'].items() if 1 < sum(k) < 7])

    axwidth = barwidth*len(bkeys)/pltwidth
    # [0.05, 0.3, 0.2226, 0.1]
    axC0 = fig2.add_axes([0.05, 0.3, 0.025, axwidth])
    axC1 = fig2.add_axes([0.085, 0.3, 0.275, axwidth])
    axC2 = fig2.add_axes([0.39, 0.3, 0.15, axwidth])
    axC3 = fig2.add_axes([0.58, 0.3, 0.15, axwidth])
    axC4 = fig2.add_axes([0.77, 0.3, 0.15, axwidth])

    crosspos = [7, 8, 13]
    checkpos = [i for i in range(len(bkeys)) if i not in crosspos]
    axC0.scatter(x=[1]*len(checkpos), y=checkpos, marker=r'$\checkmark$', color='black', s=20)
    axC0.scatter(x=[1]*len(crosspos), y=crosspos, marker=r'$\chi$', color='red', s=20)
    axC0.spines[:].set_visible(True)
    axC0.tick_params(axis='both',
                     which='both',
                     left=False,
                     top=False,
                     labelleft=False,
                     bottom=False,
                     labelbottom=False)
    axC0.set_yticks(ticks=range(15))
    axC0.set_ylim([-0.5, 14.5])
    axC0.grid(which='major', axis='y')
    axC0.set_axisbelow(True)
    axC0.set_ylabel('Matches tree branching')

    ypos, xpos = np.where(np.array(bkeys) == 1)
    xpos = 6 - xpos
    ypos = 14 - ypos
    # xpos = 6 - xpos
    axC1.scatter([x for i, x in enumerate(xpos) if ypos[i] not in crosspos], [y for i, y in enumerate(ypos) if y not in crosspos], marker='o', color='black', s=20)
    axC1.scatter([x for i, x in enumerate(xpos) if ypos[i] in crosspos], [y for i, y in enumerate(ypos) if y in crosspos], marker='o', color='black', s=20)
    axC1.set_xticks(ticks=range(7), rotation=45)
    axC1.set_yticks(ticks=range(len(bkeys)))
    axC1.set_xticklabels(labels=branchorder[::-1])
    axC1.grid()
    axC1.set_axisbelow(True)
    axC1.set_xlim([-0.5, 6.5])
    axC1.set_ylim([-0.5, 14.5])
    axC1.tick_params(axis='y',
                     which='both',
                     left=False,
                     top=False,
                     labelleft=False)
    axC1.tick_params(axis='x', rotation=45)
    axC1.set_xticks(labels=[bnamedict[i.get_text()] for i in axC1.get_xticklabels()], ticks=axC1.get_xticks(), rotation=0)
    axC1.set_xlabel('Branches with SMs')
    axC1.spines[:].set_visible(False)
    axdict = {0: axC2, 1: axC3, 2: axC4}
    for i, l in enumerate('L1 L2 shared'.split()):
        data = cnts[l]
        total_cnt = sum(data.values())
        ax = axdict[i]
        # ax.barh(y=[str(b) for b in bkeys[::-1]], width=[data[k]/total_cnt for k in bkeys[::-1]], color=colour[l], height=0.75)
        ax.barh(y=[str(b) for b in bkeys[::-1]], width=[data[k] for k in bkeys[::-1]], color=colour[l], height=0.75)
        ax.tick_params(axis='y',
                       which='both',
                       left=False,
                       top=False,
                       labelleft=False)
        ax.set_ylim([-0.5, 14.5])
        ax.set_xlabel(f'{l} SMs')
        ax.axhline(y=6.5, color='grey', linestyle='dotted')
        ax.axhline(y=13.5, color='grey', linestyle='dotted')
        ax=cleanax(ax)
        ax.grid(which='major', axis='y')
        ax.set_axisbelow(True)
        # ax.set_xlim([0, [0.25, 0.25, 0.55][i]])
        ax.set_xlim([0, 25])

    axC0.text(0.25, 1, 'E', transform=axC0.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')

    # </editor-fold>

    plt.savefig(f"fig2_v2.png", dpi=300)
    plt.savefig(f"fig2_v2.pdf", dpi=300)
    plt.close()

    # </editor-fold>


    # <editor-fold desc="Figure 3 plots">
    barwidth = 0.125
    pltwidth = 4.76
    dpi = 100
    fig3 = plt.figure(figsize=[pltwidth, 6], dpi=dpi)
    trans = mtransforms.ScaledTranslation(-5/dpi, 5/dpi, fig3.dpi_scale_trans)


    # <editor-fold desc="A: Add tree image">
    def getcntplt(ax, branch, bar_width):
        # ax = fig1.add_axes(axispos)
        y_bottom = [0]
        ymax = 25
        ax.set_xlim([-0.5, 0.5])
        ax.set_ylim([0, ymax])
        ax.set_yticks([0, ymax])
        ax.tick_params(bottom=False, labelbottom=False)
        ax.set_xlabel(bnamedict[branch], fontdict={'backgroundcolor': 'white', 'linespacing': 0, 'horizontalalignment': 'center', 'bbox': {'pad': 0, 'fc': 'white', 'ec': 'white'}})
        ax.set_yticklabels(labels=ax.get_yticklabels(), fontdict={'backgroundcolor': 'white', 'linespacing': 0, 'horizontalalignment': 'center', 'verticalalignment': 'center', 'bbox': {'pad': 0, 'fc': 'white', 'ec': 'white'}})
        ax.spines[:].set_visible(False)
        for t in ('SNV', 'Indel'):
            y = branchscnt[t][branch]
            ax.bar([0], [y], label=f'{t}', bottom=y_bottom, color=colour[('leaf', t)], width=bar_width)
            y_bottom = [y_bottom[0] + y]
        return ax
    # END

    tree = mpimg.imread(f'{cwd}/../tree_sketch_2.png')
    axA = fig3.add_axes([0.05, 0.525, 0.6, 0.425])
    axA.imshow(tree, aspect='auto')
    axA.text(0, 1, 'A', transform=axA.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    axA.spines[:].set_visible(False)
    axA.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)


    df = allsmdata.loc[allsmdata.tissue == 'leaf'].copy()
    df['type'] = 'SNV'
    df.loc[[a[0] in '+-' for a in df.alt_allele], 'type'] = 'Indel'
    branchscnt = {g[0]: defaultdict(int, g[1].branch.value_counts()) for g in df.groupby('type')}
    bar_width = 0.75
    axwidth = (barwidth*1/pltwidth)/0.6
    axA1 = getcntplt(axA.inset_axes([0.9, 0.7, axwidth, 0.1]), 'wt_1', bar_width)
    axA2 = getcntplt(axA.inset_axes([0.9, 0.25, axwidth, 0.1]), 'wt_7', bar_width)
    axA3 = getcntplt(axA.inset_axes([0.75, 0.85, axwidth, 0.1]), 'wt_18', bar_width)
    axA4 = getcntplt(axA.inset_axes([0.45, 0.85, axwidth, 0.1]), 'wt_19', bar_width)
    axA5 = getcntplt(axA.inset_axes([0.25, 0.85, axwidth, 0.1]), 'mut_15', bar_width)
    axA6 = getcntplt(axA.inset_axes([0.1, 0.85, axwidth, 0.1]), 'mut_11_1', bar_width)
    axA7 = getcntplt(axA.inset_axes([0.1, 0.5, axwidth, 0.1]), 'mut_11_2', bar_width)

    # </editor-fold>


    # # <editor-fold desc="Plot number of SMs identified in leaves">
    # df = allsmdata.loc[allsmdata.tissue == 'leaf'].copy()
    # df['type'] = 'SNV'
    # df.loc[[a[0] in '+-' for a in df.alt_allele], 'type'] = 'Indel'
    # branchscnt = {g[0]: defaultdict(int, g[1].branch.value_counts()) for g in df.groupby('type')}
    #
    # axwidth = barwidth*7/7.24
    # axB = fig3.add_axes([0.55, 0.825, axwidth, 0.125])
    # y_bottom = [0]*7
    # for t in ('SNV', 'Indel'):
    #     y = [branchscnt[t][b] for b in branches]
    #     axB.bar(branches, y, label=f'{t}', bottom=y_bottom, color=colour[('leaf', t)], width=0.75)
    #     y_bottom = [y_bottom[i] + y[i] for i in range(7)]
    # axB.set_xlabel("Branch")
    # axB.set_ylabel("SMs")
    # axB.set_yticks(range(0, 21, 5))
    # axB = cleanax(axB)
    # axB.set_xlim(-0.5, 6.5)
    # axB.tick_params(axis='x', rotation=45)
    # axB.text(0, 1, 'B', transform=axB.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    # # </editor-fold>


    # <editor-fold desc="B: Plot Number of organ specific SMs">
    legcol = {'SNV': 'black', 'Indel': 'grey'}
    df = allsmdata.copy()
    df['type'] = ['SNV' if row[2][0] not in '+-' else 'Indel' for row in df.itertuples(index=False)]
    cnts = df \
        .drop_duplicates(subset='chromosome position alt_allele'.split()) \
        .organ.value_counts().to_dict()

    axB = fig3.add_axes([0.7, 0.825, 0.3, 0.15])
    venn2(subsets=(cnts['fruit'], cnts['leaf'], cnts['shared']),
          set_labels='fruit leaf'.split(),
          set_colors=(colour['theme2'], colour['leaf']), alpha=1, ax=axB)
    # axB.text(0, 1, 'B', transform=axB.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    axB.text(0, 1, 'B', transform=axB.transAxes, fontsize=BIGGER_SIZE, va='bottom', weight='bold')

    # </editor-fold>


    # <editor-fold desc="C: Venn diagram showing organ specificity of SMs">
    def tobed(df):
        df = df['chromosome position'.split()].copy()
        df.columns = 'chromosome end'.split()
        df['start'] = df.end - 1
        df = df['chromosome start end'.split()]
        return bt.BedTool.from_dataframe(df)

    leafsms = allsmdata.loc[allsmdata.tissue == 'leaf'].drop_duplicates(subset='chromosome position alt_allele'.split()).copy()
    leafsms = tobed(leafsms)
    l1sms = allsmdata.loc[allsmdata.tissue == 'L1'].drop_duplicates(subset='chromosome position alt_allele'.split())
    l1sms = tobed(l1sms)
    l2sms = allsmdata.loc[allsmdata.tissue == 'L2'].drop_duplicates(subset='chromosome position alt_allele'.split())
    l2sms = tobed(l2sms)
    leafcnt = len(leafsms.subtract(l1sms).subtract(l2sms))
    l1cnt = len(l1sms.subtract(leafsms).subtract(l2sms))
    l2cnt = len(l2sms.subtract(leafsms).subtract(l1sms))
    leafl2cnt = len(leafsms.intersect(l2sms, u=True).subtract(l1sms))
    leafl1cnt = len(leafsms.intersect(l1sms, u=True).subtract(l2sms))
    l1l2cnt = len(l1sms.intersect(l2sms, u=True).subtract(leafsms))
    l1l2leafcnt = len(l1sms.intersect(l2sms, u=True).intersect(leafsms, u=True))
    axC = fig3.add_axes([0.7, 0.525, 0.3, 0.25])
    # venn3(subsets=(leafcnt, l2cnt, leafl2cnt, l1cnt, leafl1cnt, l1l2cnt, l1l2leafcnt),
    #       set_labels='leaf L2 L1'.split(),
    #       set_colors=(colour['leaf'], colour['L2'], colour['L1']), alpha=1, ax=axC)
    # venn3(subsets=(l1cnt, leafcnt, leafl1cnt, l2cnt, l1l2cnt, leafl2cnt, l1l2leafcnt),
    #       set_labels='Fruit:L1 leaf Fruit:L2'.split(),
    #       set_colors=(colour['L1'], colour['leaf'], colour['L2']), alpha=1, ax=axC)
    venn3(subsets=(l2cnt, leafcnt, leafl2cnt, l1cnt, l1l2cnt, leafl1cnt, l1l2leafcnt),
          set_labels='fruit:L2 leaf fruit:L1'.split(),
          set_colors=(colour['L2'], colour['leaf'], colour['L1']), alpha=1, ax=axC)
    # axC.text(0, 1, 'C', transform=axC.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')

    # </editor-fold>


    # <editor-fold desc="D: Plot genotyping calls">
    # Originally was separating out SNVs and indels, but since there is not much
    # difference between them now I draw a single plot showing them together

    # Read SMs in L1, L2 and leaves
    df = pd.read_table(allsmrc, index_col=[0, 1, 2])
    gtvals = [1, 2, 5, 10, 15, 20]
    l1comb = {'L1': '10',
              'leaf': '01',
              'L1+leaf': '11',
              'None': '00'
              }
    l2comb = {'L2': '10',
              'leaf': '01',
              'L2+leaf': '11',
              'None': '00'
              }
    lcombdict = {0: l1comb, 1: l2comb}
    lcnts = pd.DataFrame()

    for i, l in enumerate('L1 L2'.split()):
        gtcnts = pd.DataFrame()
        for gt in gtvals:
            # tissues = [l, 'leaf']
            tissues = [l]
            # dfmask = df[[f"{b}{l}" for b in branches for l in ["_l1", "_l2", '']]].copy()
            dfmask = df[[f"{b}{t}" for b in branches for t in [f"_{l.lower()}", '']]].copy()
            dfmask[dfmask < gt] = -1
            dfmask[dfmask != -1] = 1
            dfmask[dfmask == -1] = 0
            # gtdict = {'SNP': defaultdict(int), 'Indel': defaultdict(int)}
            gtdict = defaultdict(int)
            for b in branches:
                bsm = allsmdata.loc[(allsmdata.branch == b) & (allsmdata.tissue.isin(tissues))].copy()
                bsm = bsm.loc[[t[0] not in '+-' for t in bsm.alt_allele]]
                print(bsm.shape)
                bsm = bsm.set_index('chromosome position alt_allele'.split())
                # for j, t in enumerate('SNP Indel'.split()):
                # Plots SNP heatmaps
                # bdf = dfmask.loc[df.var_type == t][[f"{b}{l}" for l in ["_l1", "_l2", '']]].copy()
                # bdf = dfmask[[f"{b}{l}" for l in ["_l1", "_l2", '']]].copy()
                bdf = dfmask[[f"{b}{t}" for t in [f"_{l.lower()}", '']]].copy()
                bdf = bdf.loc[bdf.index.isin(bsm.index)]
                # bdf = bdf.loc[bdf.apply(sum, axis=1) > 0]
                bdf.sort_values(list(bdf.columns), inplace=True)
                # l = list(bdf.apply(lambda x: f'{x[0]}{x[1]}{x[2]}', axis=1))
                cmb = list(bdf.apply(lambda x: f'{x[0]}{x[1]}', axis=1))
                for k, v in lcombdict[i].items():
                    gtdict[k] += cmb.count(v)

            gtdict = pd.DataFrame({'tissue': gtdict.keys(), 'count': gtdict.values()})
            gtdict['gt'] = gt
            print(gtdict)
            # print(sum(gtdict['count']), gtdict)
            gtdict['count'] /= sum(gtdict['count'])
            gtdict['count'] *= 100
            gtcnts = pd.concat([gtcnts, gtdict])

        gtcnts['Layer'] = l
        lcnts = pd.concat([lcnts, gtcnts])
    lcnts.reset_index(inplace=True, drop=True)

    axD1 = fig3.add_axes([0.1, 0.3, 0.5, 0.175])
    pltdf = lcnts.loc[(lcnts.Layer == 'L1') & (~lcnts.tissue.isin(['None', 'leaf']))]
    axD1 = sns.lineplot(pltdf, x='gt', y='count', hue='tissue', legend=False, ax=axD1, zorder=0, palette=[colour['L1'], colour['L1+leaf']])
    axD1 = sns.scatterplot(pltdf, x='gt', y='count', style='tissue', hue='tissue', markers=['o', '^', 's', '*', '8',  'p', 'd'], ax=axD1, edgecolor='black', linewidth=0, zorder=1, palette=[colour['L1'], colour['L1+leaf']], s=20)
    axD1 = cleanax(axD1)
    plt.setp(axD1.lines, linewidth=1)
    # plt.setp(axD1.collections, sizes=[20])
    axD1.set_xticks(gtvals)
    # axD1.legend(frameon=False, ncols=1, bbox_to_anchor=[0.35, 0.6])
    axD1.legend(frameon=False, ncol=3, loc='upper center')
    axD1.set_xlabel('')
    axD1.set_ylabel('SM-branch pairs (%)')
    axD1.set_xlim(-10, 30)
    axD1.set_ylim(0, 120)
    axD1.set_facecolor('antiquewhite')
    axD1.text(0, 1, 'D', transform=axD1.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    leafcnt = len(leafsms.subtract(l1sms))
    l1cnt = len(l1sms.subtract(leafsms))
    leafl1cnt = len(leafsms.intersect(l1sms, u=True))

    leafcnt = lcnts.loc[(lcnts['gt'] == 1) & (lcnts.tissue == 'leaf') & (lcnts.Layer == 'L1'), 'count'].to_list()[0]
    l1cnt = lcnts.loc[(lcnts['gt'] == 1) & (lcnts.tissue == 'L1') & (lcnts.Layer == 'L1'), 'count'].to_list()[0]
    leafl1cnt = lcnts.loc[(lcnts['gt'] == 1) & (lcnts.tissue == 'L1+leaf') & (lcnts.Layer == 'L1'), 'count'].to_list()[0]

    axD3 = fig3.add_axes([0.1, 0.3, 0.1, 0.175])
    v = venn2(subsets=(l1cnt, leafcnt, leafl1cnt),
              set_labels='L1 leaf'.split(),
              set_colors=(colour['L1'], colour['leaf']), alpha=1, ax=axD3)
    for i in '01 10 11 A B'.split():
        v.get_label_by_id(i).set_text('')

    leafcnt = lcnts.loc[(lcnts['gt'] == 20) & (lcnts.tissue == 'leaf') & (lcnts.Layer == 'L1'), 'count'].to_list()[0]
    l1cnt = lcnts.loc[(lcnts['gt'] == 20) & (lcnts.tissue == 'L1') & (lcnts.Layer == 'L1'), 'count'].to_list()[0]
    leafl1cnt = lcnts.loc[(lcnts['gt'] == 20) & (lcnts.tissue == 'L1+leaf') & (lcnts.Layer == 'L1'), 'count'].to_list()[0]

    axD3 = fig3.add_axes([0.5, 0.3, 0.1, 0.175])
    v = venn2(subsets=(l1cnt, leafcnt, leafl1cnt),
              set_labels='L1 leaf'.split(),
              set_colors=(colour['L1'], colour['leaf']), alpha=1, ax=axD3)
    for i in '01 10 11 A B'.split():
        v.get_label_by_id(i).set_text('')

    pltdf = lcnts.loc[(lcnts.Layer == 'L2') & (~lcnts.tissue.isin(['None', 'leaf']))]
    axD2 = fig3.add_axes([0.1, 0.06, 0.5, 0.175])
    axD2 = sns.lineplot(pltdf, x='gt', y='count', hue='tissue', legend=False, ax=axD2, zorder=0, palette=[colour['L2'], colour['L2+leaf']])
    axD2 = sns.scatterplot(pltdf, x='gt', y='count', style='tissue', hue='tissue', markers=['o', '^', 's', '*', '8',  'p', 'd'], ax=axD2, edgecolor='black', linewidth=0, zorder=1, palette=[colour['L2'], colour['L2+leaf']], s=20)
    axD2 = cleanax(axD2)
    plt.setp(axD2.lines, linewidth=1)
    # plt.setp(axD2.collections, sizes=[20])
    axD2.set_xticks(gtvals)
    # axD2.legend(frameon=False, ncols=1, bbox_to_anchor=[0.35, 0.6])
    # axD2.legend(frameon=False, ncol=3)
    axD2.legend(frameon=False, ncol=3, loc='upper center')
    axD2.set_xlim(-10, 30)
    axD2.set_ylim(0, 120)
    # axD2.tick_params(axis='y', left=False, labelleft=False)
    axD2.set_xlabel('Minimum alternate reads')
    axD2.set_ylabel('SM-branch pairs (%)')
    axD2.set_facecolor('lavender')
    # axD2.text(0, 1, 'F', transform=axD2.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')

    leafcnt = lcnts.loc[(lcnts['gt'] == 1) & (lcnts.tissue == 'leaf') & (lcnts.Layer == 'L2'), 'count'].to_list()[0]
    l2cnt = lcnts.loc[(lcnts['gt'] == 1) & (lcnts.tissue == 'L2') & (lcnts.Layer == 'L2'), 'count'].to_list()[0]
    leafl1cnt = lcnts.loc[(lcnts['gt'] == 1) & (lcnts.tissue == 'L2+leaf') & (lcnts.Layer == 'L2'), 'count'].to_list()[0]

    axD4 = fig3.add_axes([0.1, 0.06, 0.1, 0.175])
    v = venn2(subsets=(l2cnt, leafcnt, leafl1cnt),
              set_labels='L2 leaf '.split(),
              set_colors=(colour['L2'], colour['leaf']), alpha=1, ax=axD4)
    for i in '01 10 11 A B'.split():
        v.get_label_by_id(i).set_text('')

    leafcnt = lcnts.loc[(lcnts['gt'] == 20) & (lcnts.tissue == 'leaf') & (lcnts.Layer == 'L2'), 'count'].to_list()[0]
    l1cnt = lcnts.loc[(lcnts['gt'] == 20) & (lcnts.tissue == 'L2') & (lcnts.Layer == 'L2'), 'count'].to_list()[0]
    leafl1cnt = lcnts.loc[(lcnts['gt'] == 20) & (lcnts.tissue == 'L2+leaf') & (lcnts.Layer == 'L2'), 'count'].to_list()[0]

    axD4 = fig3.add_axes([0.5, 0.06, 0.1, 0.175])
    v = venn2(subsets=(l2cnt, leafcnt, leafl1cnt),
              set_labels='L2 leaf '.split(),
              set_colors=(colour['L2'], colour['leaf']), alpha=1, ax=axD4)
    for i in '01 10 11 A B'.split():
        v.get_label_by_id(i).set_text('')

    # </editor-fold>

    # Test venn-diagram after genotyping
    df = pd.read_table(allsmrc, index_col=[0, 1, 2])
    dfmask = dfmask = df[[f"{b}{t}" for b in branches for t in [f"_l1", f"_l2", '']]].copy()
    dfmask[dfmask < 1] = -1
    dfmask[dfmask != -1] = 1
    dfmask[dfmask == -1] = 0
    l1sms = dfmask[[f"{b}{t}" for b in branches for t in [f"_l1"]]].copy()
    l1sms = set(l1sms.loc[l1sms.apply(sum, axis=1) > 0].index.values)
    l2sms = dfmask[[f"{b}{t}" for b in branches for t in [f"_l2"]]].copy()
    l2sms = set(l2sms.loc[l2sms.apply(sum, axis=1) > 0].index.values)
    leafsms = dfmask[[f"{b}{t}" for b in branches for t in [f""]]].copy()
    leafsms = set(leafsms.loc[leafsms.apply(sum, axis=1) > 0].index.values)
    # venn3(subsets=[l1sms, l2sms, leafsms], set_colors=[colour['L1'], colour['L2'],colour['leaf']], set_labels=['L1', 'L2', 'leaf'], alpha=1)


    # <editor-fold desc="E: Counts of Leaf/L2/shared sms per branch">
    l2leafcnt = dict()
    for b in branches:
        df = allsmdata.loc[(allsmdata.branch == b) & (allsmdata.tissue != 'L1')].copy()
        df.drop_duplicates(subset='chromosome position alt_allele'.split(), inplace=True)
        l2leafcnt[bnamedict[b]] = Counter(df.organ)
    l2leafcnt = pd.DataFrame(l2leafcnt).T
    l2leafcnt.fillna(0, inplace=True)
    l2leafcnt.reset_index(inplace=True)
    l2leafcnt = l2leafcnt.melt(id_vars='index')
    l2leafcnt.loc[l2leafcnt.variable == 'fruit', 'variable'] = 'fruit'
    l2leafcnt.sort_values(['index'], inplace=True)
    # fig = plt.figure(figsize=[5, 3])
    # ax = fig.add_subplot()
    axwidth = barwidth*21/pltwidth

    axE = fig3.add_axes([0.7, 0.06, 0.29, 0.415])
    axE = sns.barplot(data=l2leafcnt, x='value', y='index', hue='variable', hue_order=['shared', 'fruit', 'leaf'], palette=[colour['fruit'], colour['L2'], colour['leaf']], orient='h', ax=axE)
    axE = cleanax(axE)
    axE.legend(title='', frameon=False, loc='lower left', bbox_to_anchor=[0.5, 0.5])
    axE.set_xlabel('L2 somatic mutations')
    axE.set_ylabel('Branch')
    axE.set_xticks(ticks=range(0, 25, 5))
    axE.text(0, 1, 'E', transform=axE.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    # plt.tight_layout()
    # plt.savefig(f'{cwd}/../leaf_l2_shared_sm_cnt.png', dpi=300)

    # </editor-fold>


    plt.savefig(f"fig3.png", dpi=300)
    plt.savefig(f"{figout}/fig3.pdf", dpi=300)
    plt.close()
    # </editor-fold>


    # <editor-fold desc="Figure 4 plots">

    # <editor-fold desc="Read data and save summary stats">
    df = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/sahu_analysis/analysis/umap_coordinates.txt')
    df.columns = 'UMAP_1 UMAP_2 branch bc cluster'.split()
    df.loc[df.branch == 'wt1', 'branch'] = 'wt_1'
    df.loc[df.branch == 'wt19', 'branch'] = 'wt_19'

    clsanno = {'L1': (3, 10),
               'L2': (0, 1, 2, 4, 5, 6),
               'L3': (7, 8, 13, 14),
               'Other': (9, 11, 12)}
    clscols = deque()
    clscols.extend(list(zip(clsanno['L1'], sns.blend_palette([colour['L1', 'SNV'], colour['L1', 'Indel']], n_colors=len(clsanno['L1']))[::-1])))
    clscols.extend(list(zip(clsanno['L2'], sns.blend_palette([colour['L2', 'SNV'], colour['L2', 'Indel']], n_colors=len(clsanno['L2']))[::-1])))
    clscols.extend(list(zip(clsanno['L3'], sns.blend_palette([colour['leaf', 'SNV'], colour['leaf', 'Indel']], n_colors=len(clsanno['L3']))[::-1])))
    clscols.extend(list(zip(clsanno['Other'], sns.blend_palette(['lightgrey', 'grey'], n_colors=len(clsanno['Other'])))))
    clscols = dict(clscols)
    df['clscols'] = [clscols[i] for i in df.cluster]
    df['celltype'] = 'Other'
    for l in 'L1 L2 L3'.split():
        df.loc[df.cluster.isin(clsanno[l]), 'celltype'] = l

    # Save summary stats
    atlasjsondata = dict()
    atlasjsondata['Number of Total cells'] = df.shape[0]
    atlasjsondata['Number of cells in each branch'] = df.groupby('branch').bc.nunique().to_dict()
    atlasjsondata['Number of cells of each celltype'] = df.celltype.value_counts().to_dict()
    atlasjsondata['Percent of cells of each celltype'] = ((df.celltype.value_counts()/df.shape[0])*100).apply(round, ndigits=1).to_dict()
    with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/axillary_plots/atlas_sm_stats.json', 'w') as fout:
        json.dump(atlasjsondata, fout, indent=2)

    # </editor-fold>

    dpi = 100
    fig4 = plt.figure(figsize=[7.24, 8], dpi=dpi)
    trans = mtransforms.ScaledTranslation(-20/dpi, 5/dpi, fig4.dpi_scale_trans)


    # <editor-fold desc="A: Draw cell-population scatter plot">
    ax = fig4.add_axes([0.65, 0.78, 0.33, 0.2])
    for j in range(0, max(df.cluster)+1):
        bdfi = df.loc[df.cluster == j]
        if bdfi.shape[0] > 0:
            ax.text(np.median(bdfi.UMAP_1), np.median(bdfi.UMAP_2), j, weight="normal")
    ax.scatter(df.UMAP_1, df.UMAP_2, s=2, c=df.clscols, linewidth=0.1)
    ax = pretty_ticks(-10, 10, 5, ax, 'x')
    ax.set_xlabel('UMAP 1')
    ax.set_ylabel('UMAP 2')
    ax = cleanax(ax)
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.text(-0.10, 1, 'B', transform=ax.transAxes , fontsize=BIGGER_SIZE, va='bottom', weight='bold')

    axleg = fig4.add_axes([0.6, 0.64, 0.39, 0.1])
    ax = axleg
    for k, v in clsanno.items():
        for c in v:
            ax.plot([], [], linewidth=3, c=clscols[c], label=c)
    legendhandels = ax.get_legend_handles_labels()
    L1leg = ax.legend(legendhandels[0][:2], legendhandels[1][:2], title='L1/Epidermis', bbox_to_anchor=[0, 0.95], bbox_transform=ax.transAxes, loc='upper left', ncol=1, alignment='left', frameon=False, borderpad=0.2, labelspacing=0.3)
    L2leg = ax.legend(legendhandels[0][2:8], legendhandels[1][2:8], title='L2/Mesophyll', bbox_to_anchor=[0.45, 0.95], bbox_transform=ax.transAxes, loc='upper left', ncol=3, alignment='left', frameon=False, borderpad=0.2, labelspacing=0.3)
    L3leg = ax.legend(legendhandels[0][8:12], legendhandels[1][8:12], title='L3/Vascular', bbox_to_anchor=[0, 0.5], bbox_transform=ax.transAxes, loc='upper left', ncol=2, alignment='left', frameon=False, borderpad=0.2, labelspacing=0.3)
    otherleg = ax.legend(legendhandels[0][12:], legendhandels[1][12:], title='Other', bbox_to_anchor=[0.45, 0.5], bbox_transform=ax.transAxes, loc='upper left', ncol=2, alignment='left', frameon=False, borderpad=0.2, labelspacing=0.3)
    ax.add_artist(L1leg)
    ax.add_artist(L2leg)
    ax.add_artist(L3leg)
    ax.spines[:].set_visible(False)
    ax.set_xlim(0, 1)
    ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)

    # </editor-fold>

    #
    # # <editor-fold desc="OBSOLETE: C: Presence of alt reads at SM positions">
    # with open(f'{cwd}/../rna.alt_readcount.pickle', 'rb') as f:
    #     mutsfilt, rnahm, isohm = pickle.load(f).values()
    # axC = deque()
    # ax = fig4.add_axes([0.8, 0.645, 0.19, 0.33])
    # cmap = sns.blend_palette(["white", colour['theme1']], as_cmap=True)
    # ax = sns.heatmap(mutsfilt, linewidths=0.5, linecolor='lightgrey', cmap=cmap, xticklabels=False, yticklabels=mutsfilt.index, cbar_kws={'label': 'Somatic mutation', 'fraction': 0.05}, cbar=False, ax=ax)
    # # ax.hlines([8], *ax.get_xlim(), color='k')
    # # ax.vlines([7], *ax.get_ylim(), color='k')
    # ax.set_ylabel('')
    # ax.set_xlabel('')
    # axC.append(ax)
    #
    # # rnahm = getmutlist(rnapos)
    # # ax = plt.subplot2grid((4, 1), (2, 0), rowspan=1, colspan=1)
    # ax = fig4.add_axes([0.8, 0.5225, 0.19, 0.11])
    # ax = sns.heatmap(rnahm, linewidths=0.5, linecolor='lightgrey', cmap=cmap, xticklabels=False, yticklabels=rnahm.index, cbar_kws={'label': 'scRNA-Seq', 'fraction': 0.05}, cbar=False, ax=ax)
    # ax.set_ylabel('RNA-Seq')
    # ylabyc = ax.transAxes.transform([0, 0.5])[1]/100
    # ax.yaxis.set_label_coords(5.2, ylabyc, transform=fig4.dpi_scale_trans)
    # ax.set_xlabel('')
    #
    # # isohm = getmutlist(isopos)
    # # ax = plt.subplot2grid((4, 1), (3, 0), rowspan=1, colspan=1)
    # ax = fig4.add_axes([0.8, 0.4, 0.19, 0.11])
    # ax = sns.heatmap(isohm, linewidths=0.5, linecolor='lightgrey', cmap=cmap, xticklabels=False, yticklabels=isohm.index, cbar_kws={'label': 'scIso-Seq', 'fraction': 0.05}, cbar=False, ax=ax)
    # ax.set_ylabel('Iso-Seq')
    # ax.set_xlabel('SNVs')
    # ylabyc = ax.transAxes.transform([0, 0.5])[1]/100
    # ax.yaxis.set_label_coords(5.2, ylabyc, transform=fig4.dpi_scale_trans)
    # # ax.set_xticks()
    # ylabyc = axC[0].transAxes.transform([0, 1])[1]/100
    # axC[0].text(5, ylabyc, 'C', transform=fig4.dpi_scale_trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    # # </editor-fold>
    #

    # <editor-fold desc="Create summary stats">
    scjsondata = dict()
    # TODO: Fix summary statistics output
    # Get number of SNVs in four branches
    expsmpos = allsmdata.loc[[True if alt[0] not in '+-' else False for alt in allsmdata.alt_allele] & (allsmdata.branch.isin(['wt_1', 'wt_19', 'mut_11_1', 'mut_15']))]
    pos = set(zip(expsmpos.chromosome, expsmpos.position, expsmpos.alt_allele))
    scjsondata["Number of SNPs in the 4 branches"] = len(pos)

    alltrans = pd.concat([rnatrans, isotrans])
    alltrans = alltrans.loc[alltrans.position != 4335333]
    alltrans.sort_values('chromosome position alt_allele'.split(), inplace=True)
    alltrans = alltrans.merge(allsmdata, on='chromosome position alt_allele'.split(), how='left')
    alltransanno = {}
    for grp in alltrans.groupby('chromosome position alt_allele'.split()):
        alltransanno[grp[0]] = tuple(grp[1].tissue.unique())
    scjsondata["SM counts with reads"] = len(alltransanno)
    scjsondata["SM counts in tissue with reads"] = {'_'.join(k): v for k, v in Counter(alltransanno.values()).items()}

    # Count number of cells with SNVs and their layers
    bcs = pd.read_table(smbcs, header=None)
    bcs = bcs.loc[~bcs[3].isna()]
    bcs[0] = bcs[0].apply(str.lower)
    bcs.sort_values([1, 2], inplace=True)
    bcdf = pd.DataFrame()
    for i, b in enumerate(['wt_1', 'wt_19', 'mut_11_1', 'mut_15']):
        bnamebcs = bcs.loc[bcs[0] == b]
        bnamebcs = set([f'{s1}_{i+1}' for s in bnamebcs[3].to_list() for s1 in s.split(',')])
        bnamedf = df.loc[df.bc.isin(bnamebcs)]
        bcdf = pd.concat([bcdf, bnamedf])

    # Write sc seq stats
    with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/axillary_plots/sc_seq_stats.json', 'w') as fout:
        json.dump(scjsondata, fout, indent=2)

    # </editor-fold>


    # <editor-fold desc="D: Clusters with mutated cells">
    # bcs = pd.read_table(smbcs, header=None)
    # bcs[0] = bcs[0].apply(str.lower)
    # bcs.sort_values([1, 2], inplace=True)
    # bcs = bcs.loc[~bcs[3].isna()]
    smids = pd.read_excel('/local/goel/Dropbox/projects/apricot_leaf/manuscript/supplementary_tables.xlsx', sheet_name='table S4 Allele_frequency', skiprows=1)
    rsmbcs = pd.read_table(smbcs, header=None)
    rsmbcs[0] = rsmbcs[0].apply(str.lower)
    rsmbcs.sort_values([1, 2], inplace=True)
    rsmbcs = rsmbcs.loc[~rsmbcs[3].isna()]

    ismbcs = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/bcs_with_sm_reads.txt', header=None)
    ismbcs[0] = ismbcs[0].apply(str.lower)
    ismbcs.sort_values([1, 2], inplace=True)
    ismbcs = ismbcs.loc[~ismbcs[3].isna()]
    bcs = deque()
    for grp in pd.concat([rsmbcs, ismbcs]).groupby([0, 1, 2]):
        bcs.append(list(grp[0]) + [','.join(set([j for i in grp[1][3] for j in i.split(',')]))])
    bcs = pd.DataFrame(bcs)
    bcs.sort_values([1, 2])

    # bcs = bcs.loc[bcs[2].isin((11957090, 19101663))]
    scind = dict(zip(['wt_1', 'wt_19', 'mut_11_1', 'mut_15'], [1, 2, 3, 4]))
    axD = deque()
    for j, grp in enumerate(bcs.groupby([1, 2])):
        # if grp[0] not in (('CUR1G', 11957090), ('CUR6G', 19101663)):
        #     continue
        smid = smids.loc[(smids.Chromosome == grp[0][0]) & (smids.Position == grp[0][1]), 'SM ID'].tolist()[0]
        for i, b in enumerate(['wt_1', 'wt_19', 'mut_15', 'mut_11_1']):
            bdf = df.loc[df.branch == b].copy()
            bdf.clscols = [[0.75, 0.75, 0.75, 1]] * bdf.shape[0]
            bdf['zorder'] = 1
            bdf['psize'] = 0.5
            bdf['lwd'] = 0
            bcb = grp[1].loc[grp[1][0] == b, 3]

            if bcb.shape[0] > 0:
                bcb = list(grp[1].loc[grp[1][0] == b, 3])[0].replace('-1', f'-1_{scind[b]}').split(',')
                if grp[0] == ('CUR6G', 3442340) and b == 'wt_19':
                    print(grp[0], b, bcb)
                    bcb.remove('TTCATTGCAAAGACTA-1_2') # Remove the barcode which has mis-aligned reads resulting in wrong genotyping
                bdf['clscols'] = [clscols[r.cluster] if r.bc in bcb else r.clscols for r in bdf.itertuples()]
                bdf['lwd'] = [0.5 if r.bc in bcb else r.lwd for r in bdf.itertuples()]
                bdf['zorder'] = [2 if r.bc in bcb else r.zorder for r in bdf.itertuples()]
                bdf['psize'] = [30 if r.bc in bcb else 0.5 for r in bdf.itertuples()]
                bdf.sort_values('zorder', inplace=True)
            # ax = fig4.add_axes([0.06 + (i*0.235), 0.51 - (j*0.1), 0.22, 0.09])
            ax = fig4.add_axes([0.06 + (i*0.23), 0.51 - (j*0.1), 0.23, 0.1])
            ax.scatter(bdf.UMAP_1, bdf.UMAP_2, s=bdf.psize, c=bdf.clscols, linewidth=bdf.lwd, edgecolor='black')
            ax.tick_params(left=False, labelleft=False, bottom=False, labelbottom=False)
            if j == 0:
                ax.set_title(bnamedict[b].replace('B', 'Branch '), pad=8, fontdict={'fontweight': 700})
            if i == 0:
                # ax.set_ylabel(f'{grp[0][0]}:{grp[0][1]}'.replace('CUR', 'Chr').replace('G', ''))
                ax.set_ylabel(smid.replace('SM', 'SM '))
            if i != 0:
                ax.set_ylabel('')
            axD.append(ax)
            # if j == 5:
            #     ax.set_xlabel(bnamedict[b])
    for i, ax in enumerate(axD):
        ax.spines[:].set_color('grey')
        if i in (2, 4, 11, 12, 19, 20):
            # axD[i].spines[:].set_linewidth(3)
            ax.patch.set_facecolor('lightgreen')
        # else:
        # ax.patch.set_facecolor('red')
        ax.patch.set_alpha(0.2)
        ax.patch.set_zorder(0)
    axD[0].text(-0.1, 1.1, 'C', transform=axD[0].transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')



    # bdf = df.loc[df.branch == 'mut_11_1'].copy()
    # bdf['zorder'] = 0
    # bdf['s'] = 0.5
    # bdf['lwd'] = 0
    #
    # bcb = list(bcs.loc[(bcs[1] == 'CUR6G') & (bcs[2] == 19101663) & (bcs[0] == 'mut_11_1'), 3])[0].replace('-1', f'-1_3').split(',')
    # bdf['clscols'] = [r.clscols if r.bc in bcb else [0.8, 0.8, 0.8, 1]  for r in bdf.itertuples()]
    # bdf['zorder'] = [1 if r.bc in bcb else r.zorder for r in bdf.itertuples()]
    # bdf['s'] = [10 if r.bc in bcb else r.s for r in bdf.itertuples()]
    # bdf['lwd'] = [0.5 if r.bc in bcb else r.lwd for r in bdf.itertuples()]
    # bdf.sort_values('zorder', inplace=True)
    # axD = fig4.add_axes([0.75, 0.075, 0.24, 0.225])
    # axD.scatter(bdf.UMAP_1, bdf.UMAP_2, s=bdf.s, c=bdf.clscols, linewidth=bdf.lwd, edgecolor='black')
    # # ax.scatter(bdf.UMAP_1, bdf.UMAP_2, s=bdf.s, c=bdf.clscols, linewidth=bdf.lw, edgecolor='white')
    # axD = cleanax(axD)
    # axD.set_title(f'mut_11_1')
    # axD.set_xlabel('UMAP 1')
    # axD.set_ylabel('UMAP 2')
    # ylabyc = axD.transAxes.transform([0, 0.5])[1]/100
    # axD.yaxis.set_label_coords(5.2, ylabyc, transform=fig4.dpi_scale_trans)
    # ylabyc = axD.transAxes.transform([0, 1])[1]/100
    # axD.text(5, ylabyc, 'D', transform=fig4.dpi_scale_trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')

    # </editor-fold>


    plt.savefig(f"fig4_v2.png", dpi=300)
    plt.savefig(f"fig4_v2.pdf", dpi=300)
    plt.close()
    # </editor-fold>


    # <editor-fold desc="Selected supplementary plots">

    # <editor-fold desc="Distribution of SMs over the genome">
    chrsize = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.fai', header=None)
    chrsize = chrsize.iloc[:8]

    pos = deque()
    for grp in datafilter.groupby('chromosome position alt_allele'.split()):
        pos.append(list(grp[0]) + list(grp[1].type.unique()) + list(grp[1].Layer.unique()) + [grp[1].shape[0]])
    pos = pd.DataFrame(pos)
    pos.columns = 'chromosome position alt_allele type Layer nbranch'.split()
    pos = pos.loc[['CUR' in c for c in pos.chromosome]]
    # pos = datafilter.drop_duplicates(subset=['chromosome', 'position', 'alt_allele']).copy()
    # pos = pos.loc[['CUR' in c for c in pos.chromosome]]
    chroms = pd.factorize(pos.chromosome)[1]
    pos.chromosome = 7 - pd.factorize(pos.chromosome)[0]
    margin = 0.2    # Change to adjust space between L1 and L2 tracks
    pos.loc[pos.Layer == 'L1', 'chromosome'] += margin
    pos.loc[pos.Layer == 'L2', 'chromosome'] -= margin
    # pos['branch_type'] = ['wt' if 'wt' in b else 'mut' for b in pos.branch]

    fig = plt.figure(figsize=[7, 5])
    ax = fig.add_subplot()
    ax = sns.scatterplot(pos, x='position', y='chromosome', hue='type', style='nbranch', markers=['o', 'X', 'X', 'X', 'X', 's'], zorder=1, size='nbranch', alpha=0.8, palette=['#540d6e', '#ee4266'], ax=ax)
    ax.set_yticklabels([''] + list(chroms)[::-1])
    ax.set_xlim([-1000000, max(chrsize[[1]].values)[0] + 1000000])
    lab = True
    for i in range(8):
        if lab:
            for j, l in enumerate(('L1', 'L2')):
                ax.axhline(i+margin-(j*2*margin), 0, chrsize.iloc[7-i, 1]/max(chrsize[[1]].values)[0], linewidth=3, zorder=0, alpha=1, color=colour[(l, 'Indel')], label=l)
            lab = False
            continue
        for j, l in enumerate(('L1', 'L2')):
            ax.axhline(i+margin-(j*2*margin), 0, chrsize.iloc[7-i, 1]/max(chrsize[[1]].values)[0], linewidth=3, zorder=0, alpha=1, color=colour[(l, 'Indel')])
    ax.legend(frameon=False)
    ax.spines[:].set_visible(False)
    ax.tick_params(left=False)
    ax.set_ylabel(None)
    ax.set_xlabel('Position (in 10Mb)', fontdict={'fontsize':10})
    ax.set_ylabel('Chromosome', fontdict={'fontsize':10})
    ax.set_title("Distribution of layer somatic mutations in the genome", fontdict={'fontsize':10})
    plt.tight_layout(pad=0.1)
    plt.savefig(f'sm_dist_over_genome.png', dpi=300)
    plt.savefig(f'sm_dist_over_genome.pdf', dpi=300)
    plt.close()

    # </editor-fold>


    # <editor-fold desc="Read depth of filtered BAM">
    def get_read_counts(fi, fo):
        with open(fi, 'r') as fin, open(fo, 'w') as fout:
            counts = [int(line.strip().split()[3]) for line in tqdm(fin)]
            counts = Counter(counts)
            for k,v in counts.items():
                fout.write(f'{k}\t{v}\n')
        return
    # END

    tmp_bdict = dict(zip([f'B{i}' for i in range(1, 8)], 'wt7 WT_1 wt18 WT_19 MUT_15 MUT_11_1 mut11_2'.split()))
    fins = [(f"{n}/bam_read_counts_b30_q10.bt2.txt", f"{b}_read_counts.txt") for b, n in tmp_bdict.items()]
    fins2 = deque()
    for n, b in bnamedict.items():
        for l in 'l1 l2'.split():
            fins2.append([f'./layer_samples/{n}/{n}_{l}/bam_read_counts_b30_q10.bt2.txt', f"{b}_{l}_read_counts.txt"])

    with Pool(processes=7) as pool:
        pool.starmap(get_read_counts, fins)
        pool.starmap(get_read_counts, fins2)

    fig = plt.figure(figsize=[6, 8], dpi=300)
    for i, b in enumerate(tmp_bdict):
        # Get leaf plot
        l = 'leaf'
        df = pd.read_table(f"{b}_read_counts.txt", header=None)
        df = df.loc[df[0] > 0]
        df = df.loc[df[0] <= 300]
        ax = fig.add_subplot(7, 3, (i*3)+1)
        ax.bar(df[0], df[1])
        ax = cleanax(ax)
        ax.set_xlabel('Depth')
        ax.set_ylabel('Frequency')
        ax.set_title(f'{b} {l}')
        # Get layer plot
        for j, l in enumerate('l1 l2'.split()):
            df = pd.read_table(f"{b}_{l}_read_counts.txt", header=None)
            df = df.loc[df[0] > 0]
            df = df.loc[df[0] <= 300]
            ax = fig.add_subplot(7, 3, (i*3) + 2 + j)
            ax.bar(df[0], df[1])
            ax = cleanax(ax)
            ax.set_xlabel('Depth')
            ax.set_ylabel('Frequency')
            ax.set_title(f'{b} {l}')
    plt.tight_layout()
    plt.savefig(f'{cwd}/mapping_depth_distribution_across_samples.pdf')

    # </editor-fold>


    # <editor-fold desc="Heterogeneity variation in regions with SM and background variation">
    # Added sample info in the VCF
    syrivcf = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/syri.vcf'
    syrivcfnew = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/nuc_div/syri.sample.vcf'
    with open(syrivcf, 'r') as fin:
        with open(syrivcfnew, 'w') as fout:
            changeversion = True
            for line in fin:
                if changeversion:
                    fout.write('##fileformat=VCFv4.2\n')
                    changeversion = False
                    continue
                if line[:2] == "##":
                    fout.write(line)
                elif line[:2] == "#C":
                    fout.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
                    fout.write(f'{line.strip()}\tFORMAT\tORA\tCUR\n')
                elif 'SYN' in line:
                    if '<' not in line:
                        fout.write(f'{line.strip()}\tGT\t1/1\t0/0\n')
    os.chdir('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/nuc_div/')
    Popen(f'bgzip syri.sample.vcf'.split())
    Popen(f'tabix syri.sample.vcf.gz'.split())
    Popen(f'vcftools --gzvcf syri.sample.vcf.gz  --window-pi 10000 --out cur_vs_ora'.split())
    syrisyn = deque()
    with open(syriout, 'r') as f:
        for line in f:
            line = line.strip().split()
            if line[10] == 'SYNAL':
                syrisyn.append(line[:3])
    syrisyn = pd.DataFrame(syrisyn)
    syrisyn[[1, 2]] = syrisyn[[1, 2]].astype(int)
    syrisyn.columns = 'Chromosome Start End'.split()
    synpr = pr.PyRanges(syrisyn)

    # Select SNV positions
    # smpos = datafilter.loc[datafilter.type=='SNV'].drop_duplicates(subset='chromosome position alt_allele'.split())
    # Select all positions
    smpos = datafilter.drop_duplicates(subset='chromosome position alt_allele'.split())
    smtype = dict(zip(smpos.chromosome.astype('str') + "_" + smpos.position.astype('str'), smpos.type))
    smdf = smpos['chromosome position'.split()].copy()
    smdf['end'] = smpos['position']
    smdf.columns = 'Chromosome Start End'.split()
    smpr = pr.PyRanges(smdf)
    smsyn = smpr.overlap(synpr).as_df()
    smpi = pd.DataFrame()
    for row in smsyn.itertuples(index=False):
        os.chdir('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/nuc_div/')
        Popen(f'vcftools --gzvcf syri.sample.vcf.gz --chr {row[0]} --from-bp {row[1]-1100}  --to-bp {row[2]+1100} --window-pi 100 --out {row[0]}_{row[1]}'.split())
        rb = int(np.floor(row[1]/100)) * 100
        rbins = np.arange(rb-1000, rb+1200, 100)
        rdf = pd.DataFrame(zip([row[0]]*(len(rbins)-1), rbins[:-1], rbins[1:]))
        rdf[1] += 1
        rdf.columns = 'CHROM  BIN_START  BIN_END'.split()
        pidf = pd.read_table(f'{row[0]}_{row[1]}.windowed.pi')
        rdf = rdf.merge(pidf, how='left')
        rdf.fillna(0, inplace=True)
        rdf['index'] = range(-10, 11)
        rdf['position'] = f'{row[0]}_{row[1]}'
        smpi = pd.concat([smpi, rdf])
    smpi['type'] = [smtype[p] for p in smpi.position]

    # smpi = smpi.pivot(index='position', columns='index', values='PI')
    meanPI = np.mean(pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/nuc_div/cur_vs_ora.windowed.pi').PI)
    fig = plt.figure(figsize=[4, 4], dpi=300)
    ax = fig.add_subplot(2, 1, 1)
    ax = sns.pointplot(smpi.loc[smpi.type=='SNV'], x='index', y='PI', ax=ax, errorbar='ci')
    ax.set_title('SNVs')
    plt.setp(ax.lines, linewidth=1)
    plt.setp(ax.collections, sizes=[10])
    ax = cleanax(ax)
    ax.set_ylabel('Nucleotide divresity')
    ax.set_xlabel('Genomic region bin (x100 bp)')
    ax.set_ylim([0, 0.015])
    ax.axhline(meanPI, linewidth=0.5, color='black', linestyle='dashed', zorder=0)
    ax = fig.add_subplot(2, 1, 2)
    ax = sns.pointplot(smpi.loc[smpi.type=='Indel'], x='index', y='PI', ax=ax, errorbar='ci')
    ax.set_title('Indels')
    plt.setp(ax.lines, linewidth=1)
    plt.setp(ax.collections, sizes=[10])
    ax = cleanax(ax)
    ax.set_ylabel('Nucleotide divresity')
    ax.set_xlabel('Genomic region bin (x100 bp)')
    ax.set_ylim([0, 0.015])
    ax.axhline(meanPI, linewidth=0.5, color='black', linestyle='dashed', zorder=0)
    plt.tight_layout(h_pad=5)
    plt.savefig(f'{cwd}/../all_sm_nuc_diversity.png', dpi=300)
    plt.savefig(f'{cwd}/../all_sm_nuc_diversity.pdf', dpi=300)
    plt.close()

    # # Code for plotting histogram on log scale
    # _, bins = np.histogram(np.array(list(binalt.values())) + 0.0000001, bins=100)
    # logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    # plt.hist(np.array(list(binalt.values())) + 0.0000001, bins=logbins)
    # plt.xscale('log')
    # plt.yscale('log')

    # </editor-fold>


    # <editor-fold desc="I: Dendrogram for SMs in L1 and L2">
    # hierarchical clustering between branches (this would show whether the SMs fit the morphology of the tree)
    # fig = plt.figure(figsize=[3, 5], dpi=300)
    fig1 = plt.figure()
    axwidth = barwidth*10/pltwidth
    # axB1 = fig1.add_axes([0.05, 0.275, 0.2, 0.125])
    # axB1 = fig1.add_axes([0.05, 0.275, 0.2, 0.125])
    axB1 = fig1.add_subplot(2, 1, 1)
    axB2 = fig1.add_subplot(2, 1, 2)
    axdict = {0: axB1, 1: axB2}
    for i, l in enumerate('L1 L2'.split()):
        pos = datafilter.loc[datafilter.Layer.isin([l, 'shared']), ['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele']].copy()
        # pos = datafilter.loc[datafilter.Layer == l, ['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele']].copy()
        alts = {(row[0], row[1], row[3]) for row in pos.itertuples(index=False)}
        branchpos = np.zeros([7, len(alts)], dtype=int)

        for j, branch in enumerate(branches):
            bdf = pos.loc[pos.branch == branch]
            balts = {(row[0], row[1], row[3]) for row in bdf.itertuples(index=False)}
            for k, alt in enumerate(alts):
                if alt in balts:
                    branchpos[j, k] = 1
        ## Get distances for hierarchical clustering
        AM = pdist(branchpos, metric='hamming')
        # AM = pdist(branchpos, lambda u, v : 1 - (sum((u==1) & (v==1))/len(u)))
        Z = linkage(AM, method='ward')
        ax = axdict[i]
        # dendrogram(optimal_leaf_ordering(Z, AM), link_color_func=lambda x: 'black', ax=ax, leaf_font_size=SMALL_SIZE)
        dendrogram(Z, link_color_func=lambda x: 'black', ax=ax, leaf_font_size=SMALL_SIZE)
        ax.spines[:].set_visible(False)
        ax.tick_params(left=False, labelleft=False)
        # ax.set_xticks(labels=[branches[i] for i in leaves_list(optimal_leaf_ordering(Z, AM))], ticks=ax.get_xticks(), rotation=90)
        ax.set_xticks(labels=[bnamedict[branches[i]] for i in leaves_list(Z)], ticks=ax.get_xticks(), rotation=0)
        ax.set_ylabel(l)
    # axB1.text(0, 1, 'I', transform=axB1.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    plt.savefig(f'{cwd}/../sm_branch_clustering.png', dpi=300)
    plt.savefig(f'{cwd}/../sm_branch_clustering.pdf', dpi=300)
    plt.close()
    # # Get maximum spanning tree
    # # Get adjacency matrix for minimum spanning tree
    # AM2 = 30 - squareform(AM)
    # np.fill_diagonal(AM2, 0)
    # g = ig.Graph.Weighted_Adjacency(AM2, mode="undirected")
    # g.vs['name'] = branches
    # inv_weight = [1./w for w in g.es["weight"]]
    # T = g.spanning_tree(weights=inv_weight)
    # print(T.is_connected())
    #
    # ax = fig.add_subplot(2, 2, 2)
    # ig.plot(g, target=ax, layout='davidson_harel', vertex_label=g.vs["name"], edge_width=np.log1p(g.es['weight']), vertex_size=0.2)
    # ax.set_xlabel('branch connectivity')
    # ax = fig.add_subplot(2, 2, 4)
    # ig.plot(T, target=ax, layout='davidson_harel', vertex_label=T.vs["name"], edge_width=np.log1p(T.es['weight']), vertex_size=0.2)
    # ax.set_xlabel('Maximum spanning tree')
    # plt.tight_layout(h_pad=2, w_pad=3)
    # plt.savefig(f'{cwd}/sm_branch_clustering.png')
    # plt.close()

    # </editor-fold>


    # <editor-fold desc="Example plot showing mutations having higher AF in L1 in all branches">
    pos = ('CUR1G', 34111279, '+AT')
    # cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
    # branches = ['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']
    readcnt = defaultdict(dict)
    for bname, bid in bnamedict.items():
        # for l in ('l1', 'l2', 'l3'):
        for l in ('l1', 'l2'):
            with open(f'{allsmrcdir}/{bname}.{l}.all_sm_in_all_samples.b30.q10.read_count.txt', 'r') as fin:
                for line in fin:
                    line = line.strip().split()
                    if line[0] == pos[0]:
                        if int(line[1]) == pos[1]:
                            try:
                                index = line.index(pos[2])
                            except ValueError:
                                readcnt[bid][l] = 0
                                break
                            readcnt[bid][l] = int(line[index+1])
                            break
                    else:
                        readcnt[bid][l] = 0
                        break
    pltdf = pd.DataFrame(readcnt).T
    pltdf.columns = 'L1 L2'.split()
    fig = plt.figure(figsize=[5, 2], dpi=300)
    ax = fig.add_subplot()
    ax = cleanax(ax)
    ax = pltdf.plot.bar(ax=ax, color=colour)
    ax.set_ylabel("Allele Frequency")
    ax.set_title("Layer specific somatic mutation present in all branches")
    ax.legend(frameon=False)
    ax.tick_params(axis='x', rotation=0)
    ax = cleanax(ax)
    plt.savefig(f"{cwd}/../layer_conserved_somatic_mutation_allele_frequency.png")
    plt.savefig(f"{cwd}/../layer_conserved_somatic_mutation_allele_frequency.pdf")
    plt.close()
    # </editor-fold>


    # <editor-fold desc="K: Correlation of SM count against branch length and branching events count from the primary branching">
    branchlength = {'wt_1': (0.93 + 1.25 + 1.5),
                    'wt_7': (0.93 + 1.25 + 1.0),
                    'wt_18': (0.93 + 2.25),
                    'wt_19': (3.13),
                    'mut_11_1': (0.77 + 0.92 + 1.51 + 0.1),
                    'mut_11_2': (0.77 + 0.92 + 1.51 + 0.2),
                    'mut_15': (0.77 + 0.92 + 0.08 + 1.45)}
    # branchcount = {'wt_1': 3,
    #                'wt_7': 3,
    #                'wt_18': 3,
    #                'wt_19': 2,
    #                'mut_11_1': 4,
    #                'mut_11_2': 4,
    #                'mut_15': 4}
    df = datafilter.copy()
    df.set_index('chromosome position alt_allele'.split(), inplace=True)
    # Filter positions that are in branches

    grp = datafilter.groupby('chromosome position alt_allele'.split()).nunique()
    grp = grp.loc[grp.branch == 7]
    df = df.loc[~df.index.isin(grp.index)]
    df = df.reset_index()

    # branchsmcnt = dict()
    layersmcount = defaultdict(dict)
    for l in 'L1 L2'.split():
        layersmcount[l] = {branch: df.loc[df.branch==branch].loc[df.Layer.isin([l, 'shared'])].shape[0] for branch in branches}
    #
    # mutsinallL1 = allsmpivot.loc[:, [b+'_L1' for b in branches]]              # Get all positions that are in L1
    # allsmmatfilt = allsmpivot.loc[mutsinallL1.apply(sum, axis=1) != 7]        # Filter positions that are in all L1
    #
    # for branch in branches:
    #     bdf = allsmmatfilt.loc[:, [f'{branch}_{l}' for l in ['L1', 'L2', 'leaf']]]
    #     bdf = bdf.loc[(bdf != 0).any(axis=1)]
    #     branchsmcnt[branch] = bdf.shape[0]
    #     layersmcount['L1'][branch] = sum(allsmmatfilt[f'{branch}_L1'])
    #     layersmcount['L2'][branch] = sum(allsmmatfilt[f'{branch}_L2'])
    #     layersmcount['leaf'][branch] = sum(allsmmatfilt[f'{branch}_leaf'])

    # fig = plt.figure(figsize=[3, 3], dpi=300)
    # ax = fig.add_subplot(4, 2, 1)
    # ax = sns.regplot(x=[branchlength[b] for b in branches], y=[branchsmcnt[b] for b in branches], ax=ax, label='All', color=colour['any'])
    # r, p = sp.stats.spearmanr([branchlength[b] for b in branches], [branchsmcnt[b] for b in branches])
    # ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p), transform=ax.transAxes)
    # # ax.set_xlabel('branch length (in m)')
    # ax.set_ylabel('Total number of SM')
    # ax = cleanax(ax)
    # ax = fig.add_subplot(4, 2, 2)
    # ax = sns.regplot(x=[branchcount[b] for b in branches], y=[branchsmcnt[b] for b in branches], ax=ax, label='All', color=colour['any'])
    # r, p = sp.stats.spearmanr([branchcount[b] for b in branches], [branchsmcnt[b] for b in branches])
    # ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p), transform=ax.transAxes)
    # # ax.set_xlabel('Number of branching event')
    # ax = cleanax(ax)
    # for l in ['L1', 'L2', 'leaf']:
    # ax = fig.add_subplot(2, 1, i)
    # ax = fig.add_subplot()
    # 0.62, 0.075, 0.1, 0.1]
    fig = plt.figure(figsize=[4, 4])
    # axD = fig1.add_axes([0.8, 0.075, 0.2, 0.1])
    axD = fig.add_subplot()
    for i, l in enumerate(['L1', 'L2']):
        axD = sns.regplot(x=[branchlength[b] for b in branches], y=[layersmcount[l][b] for b in branches], ax=axD, label=l, color=colour[l], scatter_kws={'s': 5}, ci=None)
        r, p = spearmanr([branchlength[b] for b in branches], [layersmcount[l][b] for b in branches])
        axD.text(.6, .2 - (i/10), 'r={:.2f}, p={:.2g}'.format(r, p), transform=axD.transAxes, label=l, color=colour[l])
        axD.set_ylabel(f'SMs')
        axD = cleanax(axD)
        i += 1
        # ax = fig.add_subplot(4, 2, i)
        # ax = sns.regplot(x=[branchcount[b] for b in branches], y=[layersmcount[l][b] for b in branches], ax=ax, label=l, color=colour[l])
        # r, p = sp.stats.spearmanr([branchcount[b] for b in branches], [layersmcount[l][b] for b in branches])
        # ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p), transform=ax.transAxes)
        # ax = cleanax(ax)
        # i+=1
        # # ax.set_ylabel(f'Number of SM in {l}')
        # if l == 'leaf':
        #     ax.set_xlabel('Number of branching event')
    # plt.tight_layout(h_pad=2, w_pad=2)
    axD.legend(frameon=True).remove()
    axD.set_xlabel('Branch length (in m)')
    plt.tight_layout()
    # axD.text(0, 1, 'K', transform=axD.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    plt.savefig(f'{cwd}/../branch_sm_correlation.png', dpi=300)
    plt.savefig(f'{cwd}/../branch_sm_correlation.pdf', dpi=300)
    plt.close()
    # </editor-fold>


    # <editor-fold desc="E: AF plots of L2/leaf specific mutations">
    # Get SMs in each tissue type
    smslist = dict()
    df = allsmdata.set_index('chromosome position alt_allele'.split()).copy()
    for grp in df.groupby('tissue'):
        smslist[grp[0]] = set(grp[1].index.values)
    l2pos = smslist['L2'].difference(smslist['leaf'])
    leafpos = smslist['leaf'].difference(smslist['L2'])
    bothpos = smslist['L2'].intersection(smslist['leaf'])
    # Get afvalues for the selected positions
    afvalues = pd.read_table(f'{cwd}/../all_sm_in_all_samples.allele_freq.txt')
    afvalues.set_index(['Unnamed: 0', 'Unnamed: 1', 'Unnamed: 2'], inplace=True)

    poscnts = pd.DataFrame()
    for i, (k, v) in enumerate({'L2': l2pos, 'leaf': leafpos, 'both': bothpos}.items()):
        afcnts = deque()
        for row in afvalues.itertuples(index=True):
            if row.Index in v:
                # s = df.loc[row.Index].to_numpy()[0] if k == 'leaf' else df.loc[row.Index].to_numpy()[0][0]
                bnames = df.loc[df.tissue.isin(['L2', 'leaf'])].loc[row.Index].branch.unique()
                for bname in bnames:
                    afcnts.append(list(row.Index) + [bname] + getvalues(row, [row._fields.index(c) for c in [bname, f'{bname}_l1', f'{bname}_l2']]))
        afcnts = pd.DataFrame(afcnts)
        afcnts['Group'] = k
        # afcnts.reset_index(inplace=True)
        poscnts = pd.concat([poscnts, afcnts])
    poscnts.columns = ['chromosome', 'position', 'alt_allele', 'branch', 'leaf', 'L1', 'L2', 'Group']
    # Remove posiltion with manually curated alt-allele
    poscnts = poscnts.loc[poscnts.position != 20360183]
    # Generate the plot for SNPs
    pltdf = poscnts.loc[[i[0] not in '+-' for i in poscnts.alt_allele]]
    fig = plt.figure(figsize=[3, 6])
    # axE = fig3.add_axes([0.65,  0.29, 0.3, 0.175])
    axE = fig.add_subplot(2, 1, 1)
    axE = sns.scatterplot(pltdf, x='L2', y='leaf', hue='Group', ax=axE, linewidth=0, zorder=1, palette={'L2': colour['L2'], 'leaf': colour['leaf'], 'both': colour['L2+leaf']}, s=10)
    axE.plot([0, 1], [0, 1], linestyle='dotted', color='lightgrey', zorder=0, lw=1)
    axE.set_xlim(-0.05, 1.05)
    axE.set_ylim(-0.05, 1.05)
    axE = cleanax(axE)
    # axE.legend(frameon=False).remove()
    axE.set_title('SNVs')
    # axE.text(0, 1, 'E', transform=axE.transAxes + trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    # Generate the plot for Indels
    pltdf = poscnts.loc[[i[0] in '+-' for i in poscnts.alt_allele]]
    # fig, ax = plt.subplots()
    ax = fig.add_subplot(2, 1, 2)
    ax = sns.scatterplot(pltdf, x='L2', y='leaf', hue='Group', ax=ax, linewidth=0, zorder=1, palette={'L2': colour['L2'], 'leaf': colour['leaf'], 'both': colour['L2+leaf']}, s=10)
    ax.plot([0, 1], [0, 1], linestyle='dotted', color='lightgrey', zorder=0, lw=1)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax = cleanax(ax)
    ax.set_title('Indels')
    # ax.legend(frameon=False).remove()
    plt.tight_layout()
    plt.savefig(f'{cwd}/../l2_leaf_allele_freq.png', dpi=300)
    plt.savefig(f'{cwd}/../l2_leaf_allele_freq.pdf', dpi=300)
    plt.close()
    # </editor-fold>


    # <editor-fold desc="IMPORTANT DO NOT DELETE: Heatmap comparing allele frquencies of somatic variations in leaf, L1 and L2">
    # IMPORTANT DO NOT DELETE: This section creates the all_sm_in_all_samples.read_counts.txt and all_sm_in_all_samples.alelle_freq.txt files required for main figures
    ## get allsmmat regions for use with BAM-readcount
    df = allsmdata['chromosome position position'.split()].drop_duplicates()
    df.to_csv(f'{cwd}/all_sm_in_all_samples.manually_selected.cleaned.regions', index=False, header=False, sep='\t')
    ## Use read-counts as in all_sm_in_all_samples_readcounts (layer_specific_dna_analysis_all_samples.sh)
    pos = set(list(map(tuple, allsmdata['chromosome position alt_allele'.split()].to_numpy())))
    rc = dict()
    af = dict()
    ## Read read counts at the layer SM positions in leaves
    for b in ('WT_1', 'wt7', 'wt18', 'WT_19', 'MUT_11_1', 'mut11_2', 'MUT_15'):
        with open(f'{allsmrcdir}/{b}.all_sm_in_all_samples.b30.q10.read_count.txt', 'r') as fin:
            bname = sdict[b]
            rc[bname] = dict()
            af[bname] = dict()
            for line in fin:
                line = line.strip().split()
                alts = deque()
                for p in pos:
                    if (line[0], int(line[1])) == (p[0], p[1]):
                        alts.append(p[2])
                for alt in alts:
                    try:
                        rc[bname][line[0], int(line[1]), alt] = int(line[basedict[alt]])
                        af[bname][line[0], int(line[1]), alt] = round(int(line[basedict[alt]])/int(line[3]), 2)
                    except KeyError:
                        try:
                            i = line.index(alt)
                            rc[bname][line[0], int(line[1]), alt] = int(line[i+1])
                            af[bname][line[0], int(line[1]), alt] = round(int(line[i+1])/int(line[3]), 2)
                        except ValueError:
                            rc[bname][line[0], int(line[1]), alt] = 0
                            af[bname][line[0], int(line[1]), alt] = 0
    ## Read counts at the layer SM positions in layers
    for l in ('l1', 'l2'):
        for b in ('wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15'):
            with open(f'{allsmrcdir}/{b}.{l}.all_sm_in_all_samples.b30.q10.read_count.txt', 'r') as fin:
                bname = f'{b}_{l}'
                rc[bname] = dict()
                af[bname] = dict()
                for line in fin:
                    line = line.strip().split()
                    alts = deque()
                    for p in pos:
                        if (line[0], int(line[1])) == (p[0], p[1]):
                            alts.append(p[2])
                    for alt in alts:
                        try:
                            rc[bname][line[0], int(line[1]), alt] = int(line[basedict[alt]])
                            af[bname][line[0], int(line[1]), alt] = round(int(line[basedict[alt]])/int(line[3]), 2)
                        except KeyError:
                            try:
                                i = line.index(alt)
                                rc[bname][line[0], int(line[1]), alt] = int(line[i+1])
                                af[bname][line[0], int(line[1]), alt] = round(int(line[i+1])/int(line[3]), 2)
                            except ValueError:
                                rc[bname][line[0], int(line[1]), alt] = 0
                                af[bname][line[0], int(line[1]), alt] = 0
    allsmorgan = allsmdata['chromosome position alt_allele organ'.split()].drop_duplicates()
    allsmorgan.set_index('chromosome position alt_allele'.split(), inplace=True)
    rcvalues = pd.DataFrame(rc)
    afvalues = pd.DataFrame(af)
    rcvalues['var_type'] = ['Indel' if i[2][0] in '+-' else 'SNP' for i in rcvalues.index.values]
    rcvalues['Organ'] = allsmorgan.loc[rcvalues.index.values]
    rcvalues['Organ'] = allsmorgan.loc[rcvalues.index.values]
    assert all(rcvalues.index.values == afvalues.index.values)
    afvalues['var_type'] = rcvalues['var_type'].copy()
    afvalues['Organ'] = rcvalues['Organ'].copy()
    rcvalues.sort_values(['Organ', 'var_type'], inplace=True, ascending=True)
    afvalues.sort_values(['Organ', 'var_type'], inplace=True, ascending=True)

    annocol = pd.DataFrame({
        'branch': list(rcvalues.columns[:7]) + [c.rsplit('_', maxsplit=1)[0] for c in rcvalues.columns[7:21]],
        'Sample': np.repeat(['Leaf', 'L1', 'L2'], [7, 7, 7])
    })
    annorow = pd.DataFrame({
        'var_type': rcvalues.var_type,
        'Organ': rcvalues.Organ
    })
    rcvalues.to_csv(f'{cwd}/all_sm_in_all_samples.read_counts.txt', sep='\t')
    afvalues.to_csv(f'{cwd}/all_sm_in_all_samples.allele_freq.txt', sep='\t')
    annocol.to_csv(f'{cwd}/all_sm_in_all_samples.annocol.txt', sep='\t')
    annorow.to_csv(f'{cwd}/all_sm_in_all_samples.annorow.txt', sep='\t')

    ## Creat pheatmap in R (file: /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/R/presentation_plots.R)
    # subprocess.call(f'inkscape {cwd}/all_layer_somatic_variants.allele_freq.pdf -M {cwd}/all_layer_somatic_variants.allele_freq.emf', shell=True)
    # </editor-fold>


    # <editor-fold desc="B: Cell population percentage in cluster">
    cnts = df.groupby('branch').cluster.value_counts()
    cnts = pd.DataFrame(cnts)
    cnts.columns = ['count']
    cnts.reset_index(inplace=True)
    cnts = cnts.pivot(index='cluster', columns='branch')
    cnts.fillna(0, inplace=True)
    cnts = cnts.apply(lambda x: x/sum(x), axis=0)
    cnts = cnts.droplevel(level=0, axis=1)
    # cnts = cnts[['wt_1', 'wt_19', 'mut_11_1', 'mut_15']]
    border = ['mut_11_1', 'mut_15', 'wt_19', 'wt_1']
    cnts = cnts[border]
    fig = plt.figure(figsize=[5, 2])
    # axB = fig4.add_axes([0.725, 0.875, 0.26, 0.075])
    axB = fig.add_subplot()
    y_bot = np.array([0] * 4, dtype='float')
    for row in cnts.itertuples():
        axB.barh(y=border, width=row[1:], left=y_bot, color=clscols[row[0]], edgecolor='black', linewidth=0.5, label=row.Index)
        y_bot += row[1:]
    axB.legend(bbox_to_anchor=(1.01, 1.01), ncol=2)
    axB.set_xlabel('Cells (in %)')
    axB.set_ylabel('')
    axB.set_xlim([0, 1])
    axB = cleanax(axB)
    ylabyc = axB.transAxes.transform([0, 0.5])[1]/100
    axB.yaxis.set_label_coords(0.25, ylabyc, transform=fig.dpi_scale_trans)
    axB.set_yticks(labels=[bnamedict[i.get_text()] for i in axB.get_yticklabels()], ticks=axB.get_yticks())
    # ylabyc = axB.transAxes.transform([0, 1])[1]/100
    # axB.text(0.2, ylabyc, 'B', transform=fig4.dpi_scale_trans, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    # axB.text(-0.1, 1, 'B', transform=axB.transAxes, fontsize=BIGGER_SIZE, va='bottom', weight='bold')
    plt.savefig(f'{cwd}/../cell_population_in_cluster.png', dpi=300)
    plt.savefig(f'{cwd}/../cell_population_in_cluster.pdf', dpi=300)
    plt.close()
    # </editor-fold>


    # <editor-fold desc="C: Presence of alt reads at SM positions">
    def getsmafplot(transdata, mutsdf, ax):
        transdata['data'] = 'rna'
        transdata = transdata.loc[transdata.position != 4335333].copy()
        transdata.sort_values('chromosome position alt_allele'.split(), inplace=True)
        transum = transdata.groupby('chromosome position alt_allele branch'.split()).agg('sum')
        transum.reset_index(inplace=True)
        transum['af'] = transum.ac / transum.rc
        transum = transum.merge(mutsdf, on='chromosome position alt_allele'.split(), how='left')
        transum['Layer'] = ['L1' if t == 'L1' else 'L2/leaf' for t in transum.tissue]
        transum['Branch'] = transum.branch_x == transum.branch_y
        transum.loc[transum.Branch == True, 'Branch'] = 'Match'
        transum.loc[transum.Branch == False, 'Branch'] = 'Mismatch'
        ax = sns.scatterplot(x=transum.rc, y=transum.af, hue=transum['Layer'], palette={'L1': colour['L1'], 'L2/leaf': colour['L2']}, style=transum['Branch'], markers={'Match': 'o', 'Mismatch': 'X'}, ax=ax)
        return ax
    # END

    with open(f'{cwd}/../rna.alt_readcount.pickle', 'rb') as f:
        mutsfilt, rnahm, isohm, rnatrans, isotrans, mutsdf = pickle.load(f).values()

    fig = plt.figure(figsize=[2, 4])
    # axC1 = fig4.add_axes([0.725, 0.65, 0.26, 0.15])
    axC1 = fig.add_subplot(2, 1, 1)
    axC1 = getsmafplot(rnatrans, mutsdf, axC1)
    axC1 = cleanax(axC1)
    axC1.set_xlabel('')
    axC1.set_xlabel('scRNA-Seq Read Count')
    axC1.set_ylabel('Allele frequency')
    handles, labels = axC1.get_legend_handles_labels()
    axC1.legend(handles=getvalues(handles, [1,2,4,5]), labels=getvalues(labels, [1,2,4,5]), frameon=False, loc='upper right')
    # axC1.text(-0.1, 1, 'C', transform=axC1.transAxes, fontsize=BIGGER_SIZE, va='bottom', weight='bold')

    # axC2 = fig4.add_axes([0.725, 0.425, 0.26, 0.15])
    axC2 = fig.add_subplot(2, 1, 2)
    axC2 = getsmafplot(isotrans, mutsdf, axC2)
    axC2 = cleanax(axC2)
    axC2.legend().remove()
    axC2.set_xlabel('scISO-Seq Read Count')
    axC2.set_ylabel('Allele frequency')
    plt.savefig(f'{cwd}/../alt_reads_sm_positions.png', dpi=300)
    plt.savefig(f'{cwd}/../alt_reads_sm_positions.pdf', dpi=300)
    plt.close()

    # </editor-fold>


    # <editor-fold desc="Genomic width covered by scRNA-seq and scISO-seq at different read depths">
    tmpbnamedict = {'WT_1': 'B2',
                    'WT_19': 'B4',
                    'MUT_15': 'B5',
                    'MUT_11_1': 'B6'}
    # read scRNA read depth
    rnabasecov = defaultdict(dict)
    for s in 'WT_1 WT_19 MUT_11_1 MUT_15'.split():
        with open(f'{scrnadir}/{s}/{s}.merged.depth.hist') as fin:
            for line in fin:
                line = line.strip().split()
                if line[1] == '0': continue
                rnabasecov[s][int(line[1])] = int(line[0])
    # read scISO read depth
    isobasecov = defaultdict(dict)
    for s in 'WT_1 WT_19 MUT_11_1 MUT_15'.split():
        with open(f'{scisodir}/{s}/{s}.iso_seq.depth.hist') as fin:
            for line in fin:
                line = line.strip().split()
                if line[1] == '0': continue
                isobasecov[s][int(line[1])] = int(line[0])
    # Save figures
    fig = plt.figure(figsize=[5, 5], dpi=300)
    for i, s in enumerate('WT_1 WT_19 MUT_11_1 MUT_15'.split()):
        ax = fig.add_subplot(2, 2, i+1)
        xvalues = np.array(list(rnabasecov[s].keys()))
        yvalues = np.array(list(rnabasecov[s].values()))
        ax.plot(xvalues, yvalues, label='scRNA-seq')

        xvalues = np.array(list(isobasecov[s].keys()))
        yvalues = np.array(list(isobasecov[s].values()))
        ax.plot(xvalues, yvalues, label='scISO-seq')
        ax.set_xscale('log')
        # ax.set_yscale('log')
        ax.set_xlabel('Read depth')
        ax.set_ylabel('Bases')
        ax = cleanax(ax)
        ax.legend(frameon=False)
        ax.set_title(tmpbnamedict[s])
    plt.suptitle('Base count', fontsize=BIGGER_SIZE)
    plt.tight_layout()
    plt.savefig(f'{cwd}/../scrna_sciso_base_count.pdf')
    plt.savefig(f'{cwd}/../scrna_sciso_base_count.png')
    plt.close()

    fig = plt.figure(figsize=[5, 5], dpi=300)
    for i, s in enumerate('WT_1 WT_19 MUT_11_1 MUT_15'.split()):
        print(np.sum(np.array(list(rnabasecov[s].values()))), np.sum(np.array(list(isobasecov[s].values()))))
        ax = fig.add_subplot(2, 2, i+1)
        xvalues = np.array(list(rnabasecov[s].keys()))
        maxy = np.sum(np.array(list(rnabasecov[s].values())))
        yvalues = maxy - np.cumsum(np.array([0] + list(rnabasecov[s].values()))[:-1])
        ax.plot(xvalues, yvalues, label='scRNA-seq')

        xvalues = np.array(list(isobasecov[s].keys()))
        maxy = np.sum(np.array(list(isobasecov[s].values())))
        yvalues = maxy - np.cumsum(np.array([0] + list(isobasecov[s].values()))[:-1])
        ax.plot(xvalues, yvalues, label='scISO-seq')
        ax.set_xscale('log')
        # ax.set_yscale('log')
        ax.set_xlabel('Read depth')
        ax.set_ylabel('Bases')
        ax = cleanax(ax)
        ax.legend(frameon=False)
        ax.set_title(tmpbnamedict[s])
    plt.suptitle('Cumulative base count (number of bases with at least X reads)', fontsize=BIGGER_SIZE)
    plt.tight_layout()
    plt.savefig(f'{cwd}/../scrna_sciso_cumulative_base_count.pdf')
    plt.savefig(f'{cwd}/../scrna_sciso_cumulative_base_count.png')
    plt.close()


    # Count number of bases sequenced by scRNAseq and scISOseq
    # get positions sequenced with scRNA
    basecnts = defaultdict(dict)
    for s in 'WT_1 WT_19 MUT_11_1 MUT_15'.split():
        print(s)
        onlyrna = 0
        onlyiso = 0
        rnaiso = 0
        rnapos = defaultdict(deque)
        print('read rna')
        with open(f'{scrnadir}/{s}/{s}.merged.depth.txt') as fin:
            for line in fin:
                line = line.strip().split()
                if line[2] == '0': continue
                rnapos[line[0]].append(line[1])
        for k in rnapos.keys():
            rnapos[k] = set(rnapos[k])
        isopos = defaultdict(deque)
        print('read iso')
        with open(f'{scisodir}/{s}/{s}.iso_seq.depth.txt') as fin:
            for line in fin:
                line = line.strip().split()
                if line[2] == '0': continue
                isopos[line[0]].append(line[1])
        for k in isopos.keys():
            isopos[k] = set(isopos[k])
        allchr = set(list(rnapos.keys()) + list(isopos.keys()))
        for k in allchr:
            onlyrna += len(set(rnapos[k]).difference(isopos[k]))
            onlyiso += len(set(isopos[k]).difference(rnapos[k]))
            rnaiso += len(set(rnapos[k]).intersection(isopos[k]))
        basecnts[s] = [onlyrna, onlyiso, rnaiso]
        print(basecnts)
    basecnts = {'WT_1': [16182302, 10519084, 43420994], 'WT_19': [19386734, 8236102, 45881356], 'MUT_11_1': [16606297, 8979212, 37445450], 'MUT_15': [18748963, 9421568, 27997357]}
    basecnts = {'wt_1': [16182302, 10519084, 43420994], 'wt_19': [19386734, 8236102, 45881356], 'mut_11_1': [16606297, 8979212, 37445450], 'mut_15': [18748963, 9421568, 27997357]}
    basecntdf = pd.DataFrame(basecnts).T
    basecntdf.reset_index(inplace=True)
    basecntdf.columns = 'Branch scRNA-seq scISO-seq Both'.split()
    basecntdf = basecntdf.melt(id_vars='Branch')
    fig, ax = plt.subplots(figsize=[3, 3], dpi=300)
    ax = sns.barplot(basecntdf, x='Branch', y='value', hue='variable', order=['wt_1', 'wt_19', 'mut_15', 'mut_11_1'], ax=ax)
    ax.set_xticklabels([bnamedict[k] for k in ['wt_1', 'wt_19', 'mut_15', 'mut_11_1']])
    ax.set_xlabel('Branch')
    ax.set_ylabel('Number of bases (in $10^{7}$)')
    ax.set_ylabel('SMs (x$10^{-8}$/Mbp)')
    ax.legend(frameon=False)
    plt.tight_layout(pad=0.1)
    ax = cleanax(ax)
    plt.savefig(f'{cwd}/../scrna_sciso_base_sets.pdf')
    plt.savefig(f'{cwd}/../scrna_sciso_base_sets.png')
    plt.close()

    # </editor-fold>


    # <editor-fold desc="Plot BCs containing SMs for all SMs">
    smids = pd.read_excel('/local/goel/Dropbox/projects/apricot_leaf/manuscript/supplementary_tables.xlsx', sheet_name='table S4 Allele_frequency', skiprows=1)

    rsmbcs = pd.read_table(smbcs, header=None)
    rsmbcs[0] = rsmbcs[0].apply(str.lower)
    rsmbcs.sort_values([1, 2], inplace=True)
    rsmbcs = rsmbcs.loc[~rsmbcs[3].isna()]

    ismbcs = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/bcs_with_sm_reads.txt', header=None)
    ismbcs[0] = ismbcs[0].apply(str.lower)
    ismbcs.sort_values([1, 2], inplace=True)
    ismbcs = ismbcs.loc[~ismbcs[3].isna()]
    bcs = deque()
    for grp in pd.concat([rsmbcs, ismbcs]).groupby([0, 1, 2]):
        bcs.append(list(grp[0]) + [','.join(set([j for i in grp[1][3] for j in i.split(',')]))])
    bcs = pd.DataFrame(bcs)
    bcs.sort_values([1, 2])

    smbcslist = bcs
    scind = dict(zip(['wt_1', 'wt_19', 'mut_11_1', 'mut_15'], [1, 2, 3, 4]))
    selected = pd.DataFrame()
    fig = plt.figure(figsize=[6, 9], dpi=300)
    pltcnt = 0
    for j, grp in enumerate(bcs.groupby([1, 2])):
        smid = smids.loc[(smids.Chromosome == grp[0][0]) & (smids.Position == grp[0][1]), 'SM ID'].tolist()[0]
        print(smid)
        for i, b in enumerate(['wt_1', 'wt_19', 'mut_15', 'mut_11_1']):
            pltcnt += 1
            bdf = df.loc[df.branch == b].copy()
            bdf.clscols = [[0.75, 0.75, 0.75, 1]] * bdf.shape[0]
            bdf['zorder'] = 0
            bdf['psize'] = 0.5
            bdf['lwd'] = 0
            bcb = grp[1].loc[grp[1][0] == b, 3]
            if bcb.shape[0] > 0:
                bcb = list(grp[1].loc[grp[1][0] == b, 3])[0].replace('-1', f'-1_{scind[b]}').split(',')
                if grp[0] == ('CUR6G', 3442340) and b == 'wt_19':
                    bcb.remove('TTCATTGCAAAGACTA-1_2') # Remove the barcode which has mis-aligned reads resulting in wrong genotyping
                bdf['clscols'] = [clscols[r.cluster] if r.bc in bcb else r.clscols for r in bdf.itertuples()]
                bdf['lwd'] = [0.5 if r.bc in bcb else r.lwd for r in bdf.itertuples()]
                bdf['zorder'] = [1 if r.bc in bcb else r.zorder for r in bdf.itertuples()]
                bdf['psize'] = [30 if r.bc in bcb else 0.5 for r in bdf.itertuples()]
                tmp = bdf.loc[bdf.bc.isin(bcb)].copy()
                tmp['chromosome'] = grp[0][0]
                tmp['position'] = grp[0][1]
                selected = pd.concat([selected, tmp])     # To compare selected BCs distribution for SM and WT alleles
                # selected = pd.concat([selected, bdf.loc[bdf.bc.isin(bcb)]])     # To compare selected BCs distribution for SM and WT alleles
                bdf.sort_values('zorder', inplace=True)
            # ax = fig.add_axes([0.05 + (i*0.25), 0.87 - (j*0.17), 0.2, 0.11])
            ax = fig.add_subplot(6, 4, pltcnt)
            ax.scatter(bdf.UMAP_1, bdf.UMAP_2, s=bdf.psize, c=bdf.clscols, linewidth=bdf.lwd, edgecolor='black')
            ax = cleanax(ax)
            ax.set_xlabel('UMAP 1')
            # ax.set_title(f'{grp[0][0]} {grp[0][1]} {b}')
            ax.set_title(bnamedict[b], pad=8, fontdict={'fontweight': 700})
            if i == 0:
                # ax.set_ylabel(f'{grp[0][0]} {grp[0][1]}')
                ax.set_ylabel(smid.replace('SM', 'SM '))
            if i != 0:
                ax.set_yticks([])
                ax.spines['left'].set_visible(False)
    plt.tight_layout()
    plt.savefig(f'{cwd}/../all_sm_containing_bcs.with_iso.png')
    plt.savefig(f'{cwd}/../all_sm_containing_bcs.with_iso.pdf')
    plt.close()


    rsmbcs = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/scrna_clusters/bcs_with_wt_reads.txt', header=None)
    rsmbcs[0] = rsmbcs[0].apply(str.lower)
    rsmbcs.sort_values([1, 2], inplace=True)
    rsmbcs = rsmbcs.loc[~rsmbcs[3].isna()]

    ismbcs = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/bcs_with_wt_reads.txt', header=None)
    ismbcs[0] = ismbcs[0].apply(str.lower)
    ismbcs.sort_values([1, 2], inplace=True)
    ismbcs = ismbcs.loc[~ismbcs[3].isna()]
    bcs = deque()
    for grp in pd.concat([rsmbcs, ismbcs]).groupby([0, 1, 2]):
        bcs.append(list(grp[0]) + [','.join(set([j for i in grp[1][3] for j in i.split(',')]))])
    bcs = pd.DataFrame(bcs)
    wtbcslist = bcs

    # bcs = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/scrna_clusters/bcs_with_wt_reads.txt', header=None)

    selectedwt = pd.DataFrame()
    fig = plt.figure(figsize=[6, 9], dpi=300)
    pltcnt = 0
    for j, grp in enumerate(bcs.groupby([1, 2])):
        smid = smids.loc[(smids.Chromosome == grp[0][0]) & (smids.Position == grp[0][1]), 'SM ID'].tolist()[0]
        for i, b in enumerate(['wt_1', 'wt_19', 'mut_15', 'mut_11_1']):
            pltcnt += 1
            bdf = df.loc[df.branch == b].copy()
            bdf.clscols = [[0.75, 0.75, 0.75, 1]] * bdf.shape[0]
            bdf['zorder'] = 0
            bdf['psize'] = 0.5
            bdf['lwd'] = 0
            bcb = grp[1].loc[grp[1][0] == b, 3]
            if bcb.shape[0] > 0:
                bcb = list(grp[1].loc[grp[1][0] == b, 3])[0].replace('-1', f'-1_{scind[b]}').split(',')
                bdf['clscols'] = [clscols[r.cluster] if r.bc in bcb else r.clscols for r in bdf.itertuples()]
                bdf['lwd'] = [0.5 if r.bc in bcb else r.lwd for r in bdf.itertuples()]
                bdf['zorder'] = [1 if r.bc in bcb else r.zorder for r in bdf.itertuples()]
                bdf['psize'] = [30 if r.bc in bcb else 0.5 for r in bdf.itertuples()]
                bdf.sort_values('zorder', inplace=True)
                tmp = bdf.loc[bdf.bc.isin(bcb)].copy()
                tmp['chromosome'] = grp[0][0]
                tmp['position'] = grp[0][1]
                selectedwt = pd.concat([selectedwt, tmp])     # To compare selected BCs distribution for SM and WT alleles
                # selectedwt = pd.concat([selectedwt, bdf.loc[bdf.bc.isin(bcb)]])     # To compare selected BCs distribution for SM and WT alleles
                # print(grp[0], sum(bdf.zorder))
            # ax = fig.add_axes([0.05 + (i*0.25), 0.87 - (j*0.17), 0.2, 0.11])
            ax = fig.add_subplot(6, 4, pltcnt)
            ax.scatter(bdf.UMAP_1, bdf.UMAP_2, s=bdf.psize, c=bdf.clscols, linewidth=bdf.lwd, edgecolor='black')
            ax = cleanax(ax)
            ax.set_xlabel('UMAP 1')
            ax.set_title(bnamedict[b], pad=8, fontdict={'fontweight': 700})
            if i == 0:
                # ax.set_ylabel(f'{grp[0][0]} {grp[0][1]}')
                ax.set_ylabel(smid.replace('SM', 'SM '))
            if i != 0:
                ax.set_yticks([])
                ax.spines['left'].set_visible(False)
    plt.tight_layout()
    plt.savefig(f'{cwd}/../all_wt_containing_bcs.with_iso.png')
    plt.savefig(f'{cwd}/../all_wt_containing_bcs.with_iso.pdf')
    plt.close()

    # Only scRNA: fisher_exact([[1, 37], [15, 75]], alternative='less')
    # both scRNA and scISO: fisher_exact([[3, 48], [22, 113]], alternative='less') p-value:0.046
    # </editor-fold>


    # <editor-fold desc="Correlation of number of cells with SM and WT alleles">
    smcnts = selected.groupby('chromosome position'.split()).size()
    wtcnts = selectedwt.groupby('chromosome position'.split()).size()/4
    fig, ax = plt.subplots(figsize=[3, 3], dpi=300)
    ax.scatter(wtcnts.values, smcnts.values, color='black')
    ax.set_xlim([0, 35])
    ax.set_ylim([0, 35])
    plt.axline((0, 0), slope=1, linestyle='dashed', color='grey')
    ax.set_xlabel("Average number of cells with WT allele")
    ax.set_ylabel("Number of cells with SM allele")
    ax = cleanax(ax)
    plt.savefig(f'{cwd}/../sc_cell_count_sm_vs_wt.png')
    plt.savefig(f'{cwd}/../sc_cell_count_sm_vs_wt.pdf')
    plt.close()

    # </editor-fold>

    # </editor-fold>


    # <editor-fold desc="TEMPORARY: Plot single cell clusters">
    fig = plt.figure(figsize=[12, 12])
    for i, b in enumerate(['wt_1', 'wt_19', 'mut_11_1', 'mut_15']):
        bdf = df.loc[df.branch == b]
        ax = fig.add_subplot(8, 8, (i*16) + 1)
        ax.scatter(bdf.UMAP_1, bdf.UMAP_2, s=0.5, c=bdf.clscols)
        ax.set_title(b)
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        for j in range(0, max(df.cluster)+1):
            bdfi = bdf.loc[bdf.cluster == j]
            ax = fig.add_subplot(8, 8, (i*16) + 1 + j + 1)
            ax.scatter(bdfi.UMAP_1, bdfi.UMAP_2, s=0.5, c=bdfi.clscols)
            ax.set_title(f'cluster {j}')
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
    plt.tight_layout()
    plt.savefig(f'{cwd}/../TMP_sc_cell_clusters_separated.png', dpi=600)
    # </editor-fold>


    # <editor-fold desc="Statistical test for the low number of shared SMs">
    # Create a randomised distribution
    n = 205
    d = deque()
    for i in tqdm(range(10000)):
        v = deque()
        while len(v) < n:
            a = (np.random.choice([0, 1]), np.random.choice([0, 1]))
            if 1 in a:
                v.append(a)
        d.append(Counter(v)[(1, 1)])

    # </editor-fold>


    # <editor-fold desc="Get distibution of SMs SNPs in genomic doublets. Are CpGs enriched?">
    snppos = datafilter.loc[datafilter.type == 'SNP'].drop_duplicates(subset='chromosome position alt_allele'.split())
    # TODO: question: For the doublet context, it makes sense to select the haplotype in which the mutation actually happened? Relevant when there is a SNP/indel between the haplotypes adjacent to the SM SNP.
    # Count the total number of 2-mers
    twomers = sorted(map(''.join, product(*['ATGC']*2)))
    args_list = [Namespace(fasta=Namespace(name=refcur), kmer=kmer, canonical=False, loc=False) for kmer in twomers]
    with Pool(processes=4) as pool:
        kmercnts = pool.map(cntkmer, args_list)
    kmercnts = dict(zip(twomers, kmercnts))

    refseq = readfasta(refcur)
    dup = deque()
    for p in snppos.itertuples(index=False):
        dup.append(refseq[p.chromosome][p.position-1: p.position+1])
    snppos['dup_reg'] = dup
    # read syri snps and indels
    syrisnps = pd.read_table(syrisnppath, header=None)
    syriindels = pd.read_table(syriindpath, header=None)
    # Check position where the next position is heterozygous
    printdf(snppos.loc[snppos.position.isin(syrisnps[1] - 1)])
    #     chromosome  position    branch ref_allele alt_allele  read_count  allele_freq Layer type dup_reg
    # 19       CUR1G   3423674  mut_11_2          T          A        44.0        0.253    L1  SNP      TG
    # 257      CUR5G  10982419  mut_11_1          C          T        50.0        0.431    L1  SNP      CT
    # 290      CUR6G   3442340      wt_1          A          G        27.0        0.403    L2  SNP      AG
    # 399      CUR8G   9339930     wt_19          C          T        72.0        0.456    L1  SNP      CA
    # 400      CUR8G   9696861     wt_18          A          T        51.0        0.586    L2  SNP      AT
    # Doublet change in the orangered haplotype
    # CUR5G:10982419 CG > TG
    printdf(snppos.loc[(snppos.position.isin(syriindels[1] - 1)) | (snppos.position.isin(syriindels[1] - 2))])
    #     chromosome  position branch ref_allele alt_allele  read_count  allele_freq   Layer type dup_reg
    # 26       CUR1G   9846234   wt_1          A          T        24.0        0.343      L2  SNP      AT
    # 235      CUR4G  22366515   wt_7          G          A        40.0        0.284  shared  SNP      GT
    # No SM SNP adjacent to indels in orangered

    # Change the dup_reg for CUR5G:10982419 manually
    snppos.loc[(snppos.chromosome == "CUR5G") & (snppos.position == 10982419), 'dup_reg'] = 'CG'
    snppos['dup_tar'] = [s[1] for s in snppos['dup_reg']]
    snppos['dup_tar'] = snppos.alt_allele + snppos.dup_tar
    snppos['change'] = snppos.dup_reg + u'\u279D' + snppos.dup_tar

    allchange = deque()
    for i in 'ACGT':
        for j in 'ACGT':
            if i == j:
                continue
            for k in 'ACGT':
                allchange.append(''.join([i, k, u'\u279D', j, k]))

    # dupcntabs = defaultdict()
    dupcntnorm = defaultdict()
    for l in 'L1 L2 shared'.split():
        dup = snppos.loc[snppos.Layer == l, 'change']
        dupcnt = Counter(dup)
        # dupcntabs[l] = {k: dupcnt[k] for k in kmercnts}
        dupcntnorm[l] = {k: dupcnt[k]/kmercnts[k.split(u'\u279D')[0]] for k in allchange}
    # cpgkeys = set([k for v in dupcntnorm.values() for k in v.keys() if 'CG' in k])
    # dupletjson = defaultdict()
    # dupletjson['duplet counts'] = dupcntabs
    # dupletjson['duplet count normalised'] = {k: {k1: np.format_float_scientific(v1, unique=False, precision=3) for k1, v1 in v.items()} for k, v in dupcntnorm.items()}
    # dupletjson['cgsnpcnt'] = {l: sum([dupcntabs[l][c] for c in cpgkeys]) for l in 'L1 L2 shared'.split()}
    # dupletjson['othersnpcnt'] = {l: sum([v for k, v in dupcntabs[l].items() if k not in cpgkeys]) for l in 'L1 L2 shared'.split()}
    # with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/axillary_plots/duplet_sm_stats.json', 'w') as fout:
    #     json.dump(dupletjson, fout, indent=2)
    # dupkeys = sorted(allchange, key=lambda x: sum([dupcntnorm[l][x] for l in 'L1 L2 shared'.split()]), reverse=True)
    dupkeys = allchange
    fig = plt.figure(figsize=[7.5, 2], dpi=300)
    ax = fig.add_subplot()
    ybottom = np.zeros_like(dupkeys, dtype='float')
    xloc = [i*5 + j for i in range(12) for j in range(1,5)]
    for l in 'L1 L2 shared'.split():
        y = [dupcntnorm[l][k] for k in dupkeys]
        ax.bar(xloc, y, bottom=ybottom, width=0.9, label=l, color=colour[l])
        ybottom += y
    ax.set_xlabel('doublet context for somatic mutations')
    ax.set_ylabel('Mutation ratio')
    ax.tick_params(axis='x', rotation=90)
    ax.legend(frameon=False)
    ax.set_xticks(ticks=xloc, label=dupkeys)
    ax.set_xticklabels(dupkeys)
    ax = cleanax(ax)
    plt.savefig(f'{cwd}/snp_ratio_in_genomic_doublets.png')
    plt.close()

    # </editor-fold>


    # <editor-fold desc="Check whether the SMs present in all branches are homozygous">
    ## Check how many of these positions are covered by syri annotations
    pos = [grp[0] for grp in datafilter.loc[datafilter.Layer == 'L1'].groupby('chromosome position alt_allele'.split()) if grp[1].shape[0] == 7]
    annos = ['SYN', 'SYNAL', 'INV', 'INVAL', 'TRANS', 'TRANSAL', 'INVTR',  'INVTRAL', 'DUP', 'DUPAL', 'INVDP', 'INVDPAL']
    syridf = readsyriout(syriout)
    for p in pos:
        print(p)
        print(syridf.loc[(syridf[0] == p[0]) & (syridf[1] < p[1]) & (p[1] < syridf[2])])
    # </editor-fold>


    # <editor-fold desc="Other plots that are not too informative">



    # <editor-fold desc="Plot genotyping calls">

    # Read SMs in L1, L2 and leaves
    df = pd.read_table(allsmrc, index_col=[0, 1, 2])
    gtvals = [1, 2, 5, 10, 15, 20]
    lcomb = {'l1': '100',
             'l2'  : '010',
             'leaf' : '001',
             'l1_l2'  : '110',
             'l1_leaf' : '101',
             'l2_leaf'  : '011',
             'l1_l2_leaf': '111'
             }
    gtcnts = pd.DataFrame()
    for gt in gtvals:
        dfmask = df[[f"{b}{l}" for b in branches for l in ["_l1", "_l2", '']]].copy()
        dfmask[dfmask<=gt] = -1
        dfmask[dfmask!=-1] = 1
        dfmask[dfmask==-1] = 0
        gtdict = {'SNP': defaultdict(int), 'Indel': defaultdict(int)}
        fig = plt.figure(figsize=[9, 9])
        for i, b in enumerate(branches):
            bsm = allsmdata.loc[allsmdata.branch == b]
            bsm = bsm.set_index('chromosome position alt_allele'.split())
            for j, t in enumerate('SNP Indel'.split()):
                # Plots SNP heatmaps
                bdf = dfmask.loc[df.var_type == t][[f"{b}{l}" for l in ["_l1", "_l2", '']]].copy()
                bdf = bdf.loc[bdf.index.isin(bsm.index)]
                bdf = bdf.loc[bdf.apply(sum, axis=1) > 0]
                bdf.sort_values(list(bdf.columns), inplace=True)
                l = list(bdf.apply(lambda x: f'{x[0]}{x[1]}{x[2]}', axis=1))
                for k, v in lcomb.items():
                    gtdict[t][k] += l.count(v)
                ax = fig.add_subplot(2, 7, i+1+(j*7))
                ax = sns.heatmap(bdf, ax=ax, cbar=False, linewidth=1, linecolor='red')
                ax.set_yticklabels('')
                ax.set_ylabel(t)
        plt.suptitle(f'Genetype RC > {gt}')
        plt.tight_layout()
        plt.savefig(f'{cwd}/genotyping.rc_{gt}.png')
        plt.close()
        gtdict = pd.DataFrame(gtdict)
        gtdict['gt'] = gt
        gtdict.reset_index(names='tissue', inplace=True)
        gtdict['SNP'] = (gtdict['SNP']*100)/sum(gtdict['SNP'])
        gtdict['Indel'] = (gtdict['Indel']*100)/sum(gtdict['Indel'])
        gtcnts = pd.concat([gtcnts, gtdict])

    fig = plt.figure(figsize=[4, 5], dpi=300)
    ax = fig.add_subplot(2, 1, 1)
    ax = sns.pointplot(gtcnts, x='gt', y='SNP', hue='tissue', ax=ax)
    ax = cleanax(ax)
    ax.set_ylim([-5, 80])
    plt.setp(ax.lines, linewidth=1)
    plt.setp(ax.collections, sizes=[10])
    ax.legend(frameon=False, ncols=2)
    ax.set_xlabel('Minimum read count')
    ax.set_ylabel('SNP SMs (%)')
    ax.grid(zorder=0)
    ax.set_axisbelow(True)

    ax = fig.add_subplot(2, 1, 2)
    ax = sns.pointplot(gtcnts, x='gt', y='Indel', hue='tissue', ax=ax)
    ax = cleanax(ax)
    ax.set_ylim([-5, 80])
    ax.legend('', frameon=False)
    plt.setp(ax.lines, linewidth=1)
    plt.setp(ax.collections, sizes=[10])
    ax.set_xlabel('Minimum read count')
    ax.set_ylabel('Indel SMs (%)')
    ax.grid(zorder=0)
    ax.set_axisbelow(True)

    plt.tight_layout()
    plt.savefig(f'{cwd}/../genotyping_counts.png')
    plt.close()

    # </editor-fold>


    # </editor-fold>


    # <editor-fold desc="OBSOLETE: Get correlation of mappability and SM. Ideally, this should not be different than the background distribution">

    # TODO: Fix it to work with the updated SM list
    # Read smpos same way as in repeat/gene overlap section above
    genmap = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.genmap.E0_K51.map.bedgraph'
    smindex = 0
    smanno1 = deque()
    c, p = -1, -1
    with open(genmap, 'r') as fin:
        while smindex < smpos.shape[0]:
            while tuple(smpos.iloc[smindex, [0, 1]]) == (c, p):
                smindex += 1
            c, p = smpos.iloc[smindex, [0, 1]]
            for line in fin:
                line = line.strip().split()
                if line[0] != c: continue
                if int(line[1])-25 > p: continue
                if int(line[2])-25 < p: continue
                smanno1.append((c, p, round(float(line[3]), 2)))
                break
            smindex += 1
    genmap = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.genmap.E1_K51.map.bedgraph'
    smindex = 0
    smanno2 = deque()
    c, p = -1, -1
    with open(genmap, 'r') as fin:
        while smindex < smpos.shape[0]:
            while tuple(smpos.iloc[smindex, [0, 1]]) == (c, p):
                smindex += 1
            c, p = smpos.iloc[smindex, [0, 1]]
            for line in fin:
                line = line.strip().split()
                if line[0] != c: continue
                if int(line[1])-25 > p: continue
                if int(line[2])-25 < p: continue
                smanno2.append((c, p, round(float(line[3]), 2)))
                break
            smindex += 1
    smanno1 = pd.DataFrame(smanno1)
    smanno1.columns = ['chr', 'pos', 'genmap_E0']
    smanno2 = pd.DataFrame(smanno2)
    smanno2.columns = ['chr', 'pos', 'genmap_E1']
    smanno = smanno1.merge(smanno2, how='left', on=['chr', 'pos'])
    smanno = smpos.merge(smanno, how='left', on=['chr', 'pos'])
    smanno.fillna('-', inplace=True)
    smanno.to_csv("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/somatic_mutation_genome_mappability.txt", index=False, sep='\t')
    # </editor-fold>


    return
# END


def merge_all_SM():
    """
    Merge the list of SMs identified in leaf, layers, and all samples together.

    Currently, only the SM list merging, filtering and summarising sections of this
    function are useful. All of the plotting sections are obsolete or not used in
     the manuscript. These sections are updated and used in the mutation_spectra().
    """

    # <editor-fold desc="Define imports">
    import pandas as pd
    from matplotlib import pyplot as plt
    import numpy as np
    import seaborn as sns
    from collections import deque, defaultdict, Counter, OrderedDict
    import pybedtools as bt
    from scipy.cluster.hierarchy import dendrogram, linkage, optimal_leaf_ordering, leaves_list
    from hometools.hometools import readfasta, revcomp, canonical_kmers, Namespace, cntkmer, getvalues
    from hometools.plot import cleanax
    from multiprocessing import Pool
    from itertools import product
    import scipy as sp
    import json
    import igraph as ig
    from matplotlib.colors import ListedColormap, to_rgba
    from Bio import Phylo
    from io import StringIO
    # </editor-fold>


    # <editor-fold desc="Define default values">
    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
    syriout = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out'
    smafdir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_in_all_samples_readcounts/'
    branches = ['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']
    basedict = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
    sdict = dict(zip(('WT_1', 'wt7', 'wt18', 'WT_19', 'MUT_11_1', 'mut11_2', 'mut4', 'MUT_15'), ('wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_4', 'mut_15')))
    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    # </editor-fold">


    # <editor-fold desc="Support functions">
    def read_readcounts_file(f, pos):
        """
        f: Filename to read.
        pos: target position to read. `set` of (chromosome, position, alt_allele) tuples.
        """
        basedict = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
        with open(f, 'r') as fin:
            rc = dict()
            af = dict()
            for line in fin:
                line = line.strip().split()
                alts = deque()
                for p in pos:
                    if (line[0], int(line[1])) == (p[0], p[1]):
                        alts.append(p[2])
                for alt in alts:
                    try:
                        rc[line[0], int(line[1]), alt] = int(line[basedict[alt]])
                        af[line[0], int(line[1]), alt] = round(int(line[basedict[alt]])/int(line[3]), 2)
                    except KeyError:
                        try:
                            i = line.index(alt)
                            rc[line[0], int(line[1]), alt] = int(line[i+1])
                            af[line[0], int(line[1]), alt] = round(int(line[i+1])/int(line[3]), 2)
                        except ValueError:
                            rc[line[0], int(line[1]), alt] = 0
                            af[line[0], int(line[1]), alt] = 0
        return rc, af
    # END

    # </editor-fold>


    # <editor-fold desc="Get cleaned allsmmat file">
    layerfin = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.filtered.txt'
    leaffin = "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/high_cov_mutants_sorted.all_samples.selected.txt"

    # Read layer SM
    layersm = pd.read_table(layerfin)
    layersm['tissue'] = layersm['Layer']
    layersm.drop('Layer'.split(), axis=1, inplace=True)

    # Read leaf SM
    leafsm = pd.read_table(leaffin, header=None)
    leafsm.columns = 'chromosome position ref_allele alt_allele read_count allele_freq branch selected remark read_depth ploidy'.split()
    leafsm.drop('selected remark read_depth ploidy'.split(), axis=1, inplace=True)
    leafsm['tissue'] = 'leaf'
    leafsm['type'] = 'SNP'
    leafsm.loc[[a[0] in '+-' for a in leafsm.alt_allele], 'type'] = 'Indel'

    allsm = pd.concat([leafsm, layersm])
    allsm.sort_values(['chromosome', 'position'], inplace=True)

    # Change the annotation at positions incorrectly aligned manually
    allsm.loc[allsm.position == 20360183, 'type'] = 'SNP'
    allsm.loc[allsm.position == 20360183, 'alt_allele'] = 'G'

    allsm.loc[allsm.position == 19296318, 'alt_allele'] = 'A'
    allsm.loc[allsm.position == 19296318, 'ref_allele'] = 'C'
    allsm.loc[allsm.position == 19296318, 'type'] = 'SNP'
    allsm.loc[allsm.position == 19296318, 'position'] = 19296319

    # fix mismatching branch names
    allsm.loc[allsm.branch == 'MUT_11_1', 'branch'] = 'mut_11_1'
    allsm.loc[allsm.branch == 'mut11_2', 'branch'] = 'mut_11_2'
    allsm.loc[allsm.branch == 'WT_19', 'branch'] = 'wt_19'
    allsm.loc[allsm.branch == 'WT_1', 'branch'] = 'wt_1'
    allsm.loc[allsm.branch == 'MUT_15', 'branch'] = 'mut_15'
    allsm.loc[allsm.branch == 'mut4', 'branch'] = 'mut_4'
    allsm.loc[allsm.branch == 'wt18', 'branch'] = 'wt_18'
    allsm.loc[allsm.branch == 'wt7', 'branch'] = 'wt_7'

    allsm['branch_tissue'] = allsm.branch + '_' + allsm.tissue
    # Create table for the presence and absence of SM in different tissues
    tissues = ['leaf', 'L1', 'L2']
    allbranches = ['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15', 'mut_4']
    btiss = ['_'.join([b, a]) for a in tissues for b in allbranches]
    btiss.remove('mut_4_L1')
    btiss.remove('mut_4_L2')
    allsmmat = dict()          # All somatic mutation matrix
    for g in allsm.groupby('chromosome position alt_allele'.split()):
        allsmmat[g[0]] = {k: (1 if k in g[1].branch_tissue.values else 0) for k in btiss}
    allsmmat = pd.DataFrame(allsmmat).T
    allsmmat.reset_index(inplace=True)
    allsmmat.columns = 'chromosome position alt_allele'.split() + list(allsmmat.columns.values[3:])
    leafl2mis = deque()
    for row in allsmmat.itertuples(index=False):
        if any([row[i] != row[i+15] for i in range(3, 10)]):
            leafl2mis.append(1)
        else:
            leafl2mis.append(0)
    allsmmat['leafl2mis'] = leafl2mis
    # allsmmat.to_csv('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_in_all_samples.csv', sep=',', index=False, header=True)

    # Analyse manually to select positions that are FN in leaves or L2 by comparing them
    # Read the manually curated file and update the allsmat mat
    # Removed variants at CUR2G:367749 and CUR7G:16736100 as they were complex (unclear to be de novo or gene conversion)
    allsmmat = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_in_all_samples.manually_selected.csv', sep=',')
    colname = list(allsmmat.columns.values)
    for i in range(allsmmat.shape[0]):
        if allsmmat.iat[i, 26] == 'Y':        # If selected to add manual annotation
            for v in allsmmat.iat[i, 27].split(','):
                allsmmat.iat[i, colname.index(v)] = 1
    allsmmat.set_index('chromosome position alt_allele'.split(), inplace=True)
    allsmmat.drop(['leafl2mis', 'selected', 'samples added'], axis=1, inplace=True)
    allsmmat.drop('mut_4_leaf', axis=1, inplace=True)
    allsmmat = allsmmat.loc[allsmmat.values.sum(axis=1) != 0]
    allsmmat = allsmmat.reset_index()
    allsmmat = allsmmat.melt(id_vars='chromosome position alt_allele'.split())
    allsmmat = allsmmat.loc[allsmmat.value != 0]
    allsmmat[['branch', 'tissue']] = [s.rsplit('_', 1) for s in allsmmat.variable]
    allsmmat['organ'] = 'leaf'
    allsmmat.loc[allsmmat.tissue.isin('L1 L2'.split()), 'organ'] = 'fruit'
    allsmmat.drop('variable value'.split(), axis=1, inplace=True)
    fruitsms = allsmmat.loc[allsmmat.organ == 'fruit', 'chromosome position alt_allele'.split()].to_numpy()
    fruitpos = set(list(map(tuple, fruitsms)))
    leafsms = allsmmat.loc[allsmmat.organ == 'leaf', 'chromosome position alt_allele'.split()].to_numpy()
    leafpos = set(list(map(tuple, leafsms)))
    sharedpos = fruitpos.intersection(leafpos)
    smtype = ['shared' if (row[0], row[1], row[2]) in sharedpos else row[5] for row in allsmmat.itertuples(index=False)]
    allsmmat.organ = smtype
    allsmmat.sort_values('chromosome position'.split(), inplace=True)
    # allsmmat.to_csv('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_in_all_samples.manually_selected.cleaned.csv', sep='\t', index=False)

    # Analyse manually to select positions that are FN in leaves or L1 by comparing them
    # Get AFs at sm positions selected above all_sm_in_all_samples.manually_selected.cleaned.csv
    # Only positions not in L2 to be compared as for positions in L2, there presence/absence in L1
    # and leaf has already been tested
    allsmmat = pd.read_table(f'{cwd}/all_sm_in_all_samples.manually_selected.cleaned.csv', sep=',')
    noL2pos = {grp[0] for grp in allsmmat.groupby('chromosome position alt_allele'.split()) if not grp[1].tissue.eq('L2').any()}
    # select positions with alt_AF>0.1 and alt_RC>10 for manual checking
    rc = dict()
    af = dict()
    # Read counts in leaves
    for b in ('WT_1', 'wt7', 'wt18', 'WT_19', 'MUT_11_1', 'mut11_2', 'MUT_15'):
        rc[(sdict[b], 'leaf')], af[(sdict[b], 'leaf')] = read_readcounts_file(f'{smafdir}{b}.all_sm_in_all_samples.b30.q10.read_count.txt', noL2pos)
    # Read counts in L1
    for b in ('wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15'):
        rc[(b, 'l1')], af[(b, 'l1')] = read_readcounts_file(f'{smafdir}{b}.l1.all_sm_in_all_samples.b30.q10.read_count.txt', noL2pos)
    rcdf = pd.DataFrame(rc)
    afdf = pd.DataFrame(af)
    # Leaf SMs with high AF in L1
    leafpos = {grp[0] for grp in allsmmat.groupby('chromosome position alt_allele'.split()) if grp[1].tissue.eq('leaf').all()}
    leafpos = pd.Index(leafpos)
    leafpos = rcdf.loc[leafpos]
    # No leaf only position has high read count in L1 ==> No false negatives in L1 calling
    # L1 SMs with high AF in leaf
    l1pos = {grp[0] for grp in allsmmat.groupby('chromosome position alt_allele'.split()) if grp[1].tissue.eq('L1').all()}
    l1pos = pd.Index(l1pos)
    l1pos = rcdf.loc[l1pos]
    l1pos = l1pos.loc[:, pd.Index([(b, 'leaf') for b in branches])]
    l1pos = l1pos.loc[l1pos.apply(lambda x: any(x >= 10), axis=1)]
    l1pos = afdf.loc[l1pos.index]
    # l1pos2 = l1pos.loc[:, pd.Index([(b, 'leaf') for b in branches])]
    l1pos = l1pos.loc[l1pos .loc[:, pd.Index([(b, 'leaf') for b in branches])] \
                            .apply(lambda x: any(x >= 0.1), axis=1)]
    l1pos.index.set_names(names='chromosome position alt_allele'.split(), inplace=True)
    l1possamp = {g[0]: list(g[1]) for g in allsmmat.set_index(keys='chromosome position alt_allele'.split(), drop=False) \
        .loc[l1pos.index] \
        .reset_index(drop=True) \
        .groupby('chromosome position alt_allele'.split()) \
        .branch}
    l1posfilt = pd.DataFrame()
    for grp in l1pos.groupby('chromosome position alt_allele'.split()):
        s = l1possamp[grp[0]]
        # Check if for any sample AF > 0.1 in the same sample in which L1 SM is called
        if grp[1].loc[:, list(zip(s, ['leaf']*len(s)))].apply(lambda x: any(x >= 0.1), axis=1).to_numpy()[0]:
            l1posfilt = pd.concat([l1posfilt, grp[1].loc[:, s].melt(ignore_index=False)])
    l1posfilt.columns = 'branch layer alt_freq'.split()
    # Save the positions for manual curation
    l1posfilt.to_csv('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_in_all_samples.l1_in_leaf.csv', sep='\t', index=True)
    # After manual checking no position was selected as all the candidates were in microsattelite regions
    # in the scWGS data set which had some amplification bias
    # </editor-fold>


    allsmmat = pd.read_table(f'{cwd}/all_sm_in_all_samples.manually_selected.cleaned.csv', sep=',')


    # <editor-fold desc="Create summary stats">
    allsmjsondata = dict()
    allsmjsondata["Number of events"] = allsmmat.shape[0]
    allsmjsondata["Total number of unique SMs"] = allsmmat.drop_duplicates(subset='chromosome position alt_allele'.split()).shape[0]
    allsmjsondata["Number of SNPs"] = allsmmat.loc[[c[0] not in '+-' for c in allsmmat.alt_allele]].drop_duplicates(subset='chromosome position alt_allele'.split()).shape[0]
    allsmjsondata["Number of Indels"] = allsmmat.loc[[c[0] in '+-' for c in allsmmat.alt_allele]].drop_duplicates(subset='chromosome position alt_allele'.split()).shape[0]
    allsmjsondata["Number of SMs in L1"] = allsmmat.loc[allsmmat.tissue == 'L1'].drop_duplicates(subset='chromosome position alt_allele'.split()).shape[0]
    allsmjsondata["Number of SMs in L2"] = allsmmat.loc[allsmmat.tissue == 'L2'].drop_duplicates(subset='chromosome position alt_allele'.split()).shape[0]
    allsmjsondata["Number of SMs in leaf"] = allsmmat.loc[allsmmat.tissue == 'leaf'].drop_duplicates(subset='chromosome position alt_allele'.split()).shape[0]
    allsmjsondata["Number of SMs in only fruit"] = allsmmat.loc[allsmmat.organ == 'fruit'].drop_duplicates(subset='chromosome position alt_allele'.split()).shape[0]
    allsmjsondata["Number of SMs in only leaf"] = allsmmat.loc[allsmmat.organ == 'leaf'].drop_duplicates(subset='chromosome position alt_allele'.split()).shape[0]
    allsmjsondata["Number of SMs shared"] = allsmmat.loc[allsmmat.organ == 'shared'].drop_duplicates(subset='chromosome position alt_allele'.split()).shape[0]
    allsmjsondata['Number of branches for each SM'] = OrderedDict(sorted(Counter([len(set(grp[1].branch)) for grp in allsmmat.groupby(['chromosome', 'position', 'alt_allele'])]).items()))
    with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/axillary_plots/all_sm_stats.json', 'w') as fout:
        json.dump(allsmjsondata, fout, indent=2)

    # </editor-fold>


    # <editor-fold desc="Get syri plot with SMs">
    for grp in allsmmat.groupby('organ'):
        pos = grp[1].drop_duplicates(subset='chromosome position alt_allele'.split())
        pos = pos['chromosome position'.split()]
        pos.columns = 'chromosome end'.split()
        pos['start'] = pos.end - 1
        pos = pos['chromosome start end'.split()]
        pos.to_csv(f'{cwd}/{grp[0]}_sm_pos.bed', index=False, header=False, sep='\t')

    ## Run plotsr using the above BED files as tracks (/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run)
    ## plotsr command:
    ## plotsr --sr syri.out --genomes genomes.txt --tracks tracks.txt -o syri_sm_tracks.png -S 0.3 -R -W 6 -H 8 -s 25000 -f 10
    # </editor-fold>


    # <editor-fold desc="Get syntenic and SR regions from syri.out for sm_overlap_with_structural_annotation">
    curanno = deque()
    with open(syriout, 'r') as f:
        for line in f:
            line = line.strip().split()
            if line[10] in {'SYN', 'INV', 'TRANS', 'INVTR', 'DUP', 'INVDP'}:
                curanno.append(line[:3] + [line[10]])
                continue
            if line[10] == 'NOTAL':
                if line[0] != '-':
                    curanno.append(line[:3] + [line[10]])
                # break
    curanno = pd.DataFrame(curanno)
    curanno.columns = 'chromosome start end anno'.split()
    curanno.loc[:, 'start end'.split()] = curanno.loc[:, 'start end'.split()].astype(int)
    cursyn = curanno.loc[curanno['anno'] == 'SYN']
    cursyn = bt.BedTool.from_dataframe(cursyn)
    cursr = curanno.loc[curanno['anno'] != 'SYN']
    cursr = bt.BedTool.from_dataframe(cursr)

    cursynlen = 0
    cursrlen = 0
    for grp in curanno.groupby('chromosome'):
        cursynlen += sumranges(mergeRanges(np.array(grp[1].loc[grp[1]['anno'] == 'SYN']['start end'.split()])))
        cursrlen += sumranges(mergeRanges(np.array(grp[1].loc[grp[1]['anno'] != 'SYN']['start end'.split()])))

    smdist = deque()
    for organ in 'fruit leaf shared'.split():
        sms = bt.BedTool(f'{cwd}/{organ}_sm_pos.bed').saveas()
        sms = sms.filter(func=lambda x: 'cur' not in x[0]).saveas()
        synpos = sms.intersect(cursyn, u=True)
        srpos = sms.intersect(cursr, u=True)
        smdist.append([organ, round(len(synpos)/len(sms), 2), round(len(srpos)/len(sms), 2), np.format_float_scientific((len(synpos)/len(sms))/cursynlen, unique=False, precision=3), np.format_float_scientific((len(srpos)/len(sms))/cursrlen, unique=False, precision=3)])

    smdist = pd.DataFrame(smdist)
    smdist = smdist[[0, 3, 4]]
    smdist.columns = 'organ syntenic rearranged'.split()
    smdist = pd.melt(smdist, id_vars='organ')
    smdist.value = smdist.value.astype('float')

    fig = plt.figure(figsize=[2, 2], dpi=300)
    ax = fig.add_subplot()
    ax = sns.barplot(data=smdist, x='variable', y='value', hue='organ', hue_order='fruit leaf shared'.split(), palette=colour, ax=ax)
    ax.set_ylim([0, 1.e-8])
    ax.set_xlabel('')
    ax.set_ylabel('Ratio of SMs')
    ax = cleanax(ax)
    ax.legend(frameon=False, title='')
    plt.tight_layout(pad=0.1)
    plt.savefig(f'{cwd}/sm_overlap_with_structural_annotation.png')
    plt.close()
    # </editor-fold>



    # TODO: Fix this plot so that it works for data from all_sm
    # afvalues.columns = ['chromosome', 'position', 'alt_allele'] + list(afvalues.columns)[3:]
    # fig = plt.figure(figsize=[8, 12])
    # plt.rc('font', size=10)
    # for i, grp in enumerate(datafilter.groupby('branch')):
    #     bmut = grp[1].merge(afvalues, on=['chromosome', 'position', 'alt_allele'])
    #     assert(bmut.shape[0] == grp[1].shape[0])
    #     print(grp[0], grp[1].shape[0])
    #     pltdf = bmut[['chromosome', 'position', 'Layer', 'type'] + [grp[0]+s for s in ['', '_l1', '_l2', '_l3']]]
    #     pltdf.columns = ['chromosome', 'position', 'Layer', 'type'] + ['leaf', 'l1', 'l2', 'l3']
    #     pltdf = pltdf.melt(id_vars=['chromosome', 'position', 'Layer', 'type'])
    #     ax = fig.add_subplot(4, 2, i+1)
    #     sns.lineplot(pltdf, x='variable', y='value', units='position', estimator=None, hue='Layer', style='type', ax=ax)
    #     ax.spines['right'].set_visible(False)
    #     ax.spines['top'].set_visible(False)
    #     ax.set_xlabel(grp[0])
    #     ax.set_ylabel('Alelle freq')
    #     ax.set_ylim([-0.1, 1.1])
    #     ax.get_legend().set_visible(False)
    # ax.legend(bbox_to_anchor=(1.01, 1))
    # plt.suptitle("Allele frequency of layer specific mutations")
    # plt.tight_layout(rect=(0, 0, 1, 0.98))
    # plt.savefig(f'{cwd}/af_dist_layer_sm_in_leaf.png')
    # plt.close()


    # <editor-fold desc="Get dendrogram of branches based on SM clustering">
    # Get hierarchical clustering between branches (this would show whether the
    # SMs fit the morphology of the tree). Do separately for leaf, fruit and all
    fig = plt.figure(figsize=(2, 5), dpi=300)
    for cnt, l in enumerate('leaf fruit'.split()):
        ax = fig.add_subplot(3, 1, cnt+1)
        df = allsmmat.loc[allsmmat.organ.isin([l, 'both'])]
        AM = deque()
        for i, s in enumerate(branches):
            for j, e in enumerate(branches):
                if i >= j:
                    continue
                posa = set(list(map(tuple, df.loc[df.branch == s, 'chromosome position alt_allele'.split()].to_numpy())))
                posb = set(list(map(tuple, df.loc[df.branch == e, 'chromosome position alt_allele'.split()].to_numpy())))
                # diff = sum(allsmmat[f'{s}_{l}'] != allsmmat[f'{e}_{l}'])
                diff = len(posa ^ posb)         # length of positions that are different between posa and posb (^ == XOR)
                AM.append(diff)
        AM = np.array(list(AM))
        Z = linkage(AM, method='single')
        dendrogram(optimal_leaf_ordering(Z, AM), link_color_func=lambda x: 'black', ax=ax, leaf_font_size=SMALL_SIZE)
        ax.spines[:].set_visible(False)
        ax.tick_params(left=False, labelleft=False)
        ax.set_xticks(labels=[branches[i] for i in leaves_list(optimal_leaf_ordering(Z, AM))], ticks=plt.xticks()[0], rotation=90)
        ax.set_title(l)
    # Combined
    ax = fig.add_subplot(3, 1, 3)
    df = allsmmat.copy()
    AM = deque()
    for i, s in enumerate(branches):
        for j, e in enumerate(branches):
            if i >= j:
                continue
            posa = set(list(map(tuple, df.loc[df.branch == s, 'chromosome position alt_allele'.split()].to_numpy())))
            posb = set(list(map(tuple, df.loc[df.branch == e, 'chromosome position alt_allele'.split()].to_numpy())))
            diff = len(posa ^ posb)         # length of positions that are different between posa and posb (^ == XOR)
            AM.append(diff)
    AM = np.array(list(AM))
    Z = linkage(AM, method='single')
    dendrogram(optimal_leaf_ordering(Z, AM), ax=ax, link_color_func=lambda x: 'black', leaf_font_size=SMALL_SIZE)
    ax.spines[:].set_visible(False)
    ax.tick_params(left=False, labelleft=False)
    ax.set_xticks(labels=[branches[i] for i in leaves_list(optimal_leaf_ordering(Z, AM))], ticks=plt.xticks()[0], rotation=90)
    ax.set_title('Combined')
    plt.tight_layout(h_pad=1, w_pad=1, pad=0.1)
    plt.savefig(f'{cwd}/dendrogram_for_tree_structure.png')
    plt.close()
    # </editor-fold>


    # <editor-fold desc="Get distibution of SMs SNPs in genomic triplets. Are CpGs enriched?">
    refcur = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
    df = allsmmat.copy()
    df['type'] = ['Indel' if v[0] in '+-' else 'SNP' for v in df.alt_allele]
    snppos = df.loc[df.type == 'SNP'].drop_duplicates(subset='chromosome position alt_allele'.split())
    # Count the total number of kmers
    threemers = sorted(canonical_kmers(3))
    args_list = [Namespace(fasta=Namespace(name=refcur), kmer=kmer, canonical=True) for kmer in threemers]
    with Pool(processes=4) as pool:
        kmercnts = pool.map(cntkmer, args_list)
    kmercnts = dict(zip(threemers, kmercnts))

    refseq = readfasta(refcur)
    trip = deque()
    for p in snppos.itertuples(index=False):
        trip.append(refseq[p.chromosome][p.position-2: p.position+1])
    snppos['trip_reg'] = trip
    tripcntnorm = defaultdict()
    for l in 'leaf fruit shared'.split():
        trip = snppos.loc[snppos.organ == l, 'trip_reg']
        tripcnt = Counter(trip)
        tripkeys = list(tripcnt.keys())
        tripcntclean = defaultdict(int)
        for k in tripkeys:
            kr = revcomp(k)
            kmin = min(k, kr)
            if kmin in tripcntclean:
                continue
            tripcntclean[kmin] = 0
            try:
                tripcntclean[kmin] += tripcnt[k]
                tripcnt.pop(k)
            except KeyError:
                pass
            try:
                tripcntclean[kmin] += tripcnt[kr]
                tripcnt.pop(kr)
            except KeyError:
                pass
        tripcntnorm[l] = {k: tripcntclean[k]/kmercnts[k] for k in kmercnts}
    tripkeys = sorted(kmercnts, key=lambda k: sum([tripcntnorm[l][k] for l in 'leaf fruit shared'.split()]), reverse=True)
    fig = plt.figure(figsize=[6, 3], dpi=300)
    ax = fig.add_subplot()
    ybottom = np.zeros_like(tripkeys, dtype='float')
    for l in 'leaf fruit shared'.split():
        y = [tripcntnorm[l][k] for k in tripkeys]
        ax.bar(tripkeys, y, bottom=ybottom, width=0.9, label=l, color=colour[l])
        ybottom += y
    ax.set_xlabel('Triplet context for somatic mutations')
    ax.set_ylabel('Mutation ratio')
    ax.tick_params(axis='x', rotation=90)
    ax.legend(frameon=False)
    ax = cleanax(ax)
    plt.savefig(f'{cwd}/snp_ratio_in_genomic_triplets.png')
    plt.close()
    # </editor-fold >


    # <editor-fold desc="Correlation of SM count against branch length and branching events count from the primary branching">
    branchlength = {'wt_1': (0.93 + 1.25 + 1.5),
                    'wt_7': (0.93 + 1.25 + 1.0),
                    'wt_18': (0.93 + 2.25),
                    'wt_19': (3.13),
                    'mut_11_1': (0.77 + 0.92 + 1.51 + 0.1),
                    'mut_11_2': (0.77 + 0.92 + 1.51 + 0.2),
                    'mut_15': (0.77 + 0.92 + 0.08 + 1.45)}
    branchcount = {'wt_1': 3,
                   'wt_7': 3,
                   'wt_18': 3,
                   'wt_19': 2,
                   'mut_11_1': 4,
                   'mut_11_2': 4,
                   'mut_15': 4}
    branchsmcnt = dict()
    layersmcount = defaultdict(dict)

    mutsinallL1 = allsmpivot.loc[:, [b+'_L1' for b in branches]]              # Get all positions that are in L1
    allsmmatfilt = allsmpivot.loc[mutsinallL1.apply(sum, axis=1) != 7]        # Filter positions that are in all L1

    for branch in branches:
        bdf = allsmmatfilt.loc[:, [f'{branch}_{l}' for l in ['L1', 'L2', 'leaf']]]
        bdf = bdf.loc[(bdf != 0).any(axis=1)]
        branchsmcnt[branch] = bdf.shape[0]
        layersmcount['L1'][branch] = sum(allsmmatfilt[f'{branch}_L1'])
        layersmcount['L2'][branch] = sum(allsmmatfilt[f'{branch}_L2'])
        layersmcount['leaf'][branch] = sum(allsmmatfilt[f'{branch}_leaf'])

    fig = plt.figure(figsize=[6, 10], dpi=300)
    ax = fig.add_subplot(4, 2, 1)
    ax = sns.regplot(x=[branchlength[b] for b in branches], y=[branchsmcnt[b] for b in branches], ax=ax, label='All', color=colour['any'])
    r, p = sp.stats.spearmanr([branchlength[b] for b in branches], [branchsmcnt[b] for b in branches])
    ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p), transform=ax.transAxes)
    # ax.set_xlabel('branch length (in m)')
    ax.set_ylabel('Total number of SM')
    ax = cleanax(ax)
    ax = fig.add_subplot(4, 2, 2)
    ax = sns.regplot(x=[branchcount[b] for b in branches], y=[branchsmcnt[b] for b in branches], ax=ax, label='All', color=colour['any'])
    r, p = sp.stats.spearmanr([branchcount[b] for b in branches], [branchsmcnt[b] for b in branches])
    ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p), transform=ax.transAxes)
    # ax.set_xlabel('Number of branching event')
    ax = cleanax(ax)
    i=3
    for l in ['L1', 'L2', 'leaf']:
        ax = fig.add_subplot(4, 2, i)
        ax = sns.regplot(x=[branchlength[b] for b in branches], y=[layersmcount[l][b] for b in branches], ax=ax, label=l, color=colour[l])
        r, p = sp.stats.spearmanr([branchlength[b] for b in branches], [layersmcount[l][b] for b in branches])
        ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p), transform=ax.transAxes)
        ax.set_ylabel(f'Number of SM in {l}')
        if l == 'leaf':
            ax.set_xlabel('branch length (in m)')
        ax = cleanax(ax)
        i+=1
        ax = fig.add_subplot(4, 2, i)
        ax = sns.regplot(x=[branchcount[b] for b in branches], y=[layersmcount[l][b] for b in branches], ax=ax, label=l, color=colour[l])
        r, p = sp.stats.spearmanr([branchcount[b] for b in branches], [layersmcount[l][b] for b in branches])
        ax.text(.05, .8, 'r={:.2f}, p={:.2g}'.format(r, p), transform=ax.transAxes)
        ax = cleanax(ax)
        i+=1
        # ax.set_ylabel(f'Number of SM in {l}')
        if l == 'leaf':
            ax.set_xlabel('Number of branching event')
    plt.tight_layout(h_pad=2, w_pad=2)
    plt.savefig(f'{cwd}/branching_stat_vs_sm_count.png')
    plt.close()
    # </editor-fold>


    # <editor-fold desc="Check the relationships between branching/branch lengths to SMs">
    df = allsmmat.copy()
    df.branch = pd.Categorical(df.branch)
    hassm = dict()
    for grp in df.groupby(['chromosome', 'position', 'alt_allele']):
        hassm[grp[0]] = {b: (1 if b in set(grp[1].branch) else 0) for b in branches}
    hassm = pd.DataFrame(hassm).T.reset_index(level=[0, 1, 2])
    hassm.fillna(0, inplace=True)

    branchingjson = defaultdict()
    branchingjson['mut11.1 vs mut11.2'] = sum(abs(hassm.mut_11_1 - hassm.mut_11_2))
    branchingjson['(mut11.1 & mut11.2) vs mut15'] = sum(abs((hassm.mut_11_1 * hassm.mut_11_2) - hassm.mut_15))
    branchingjson['wt1 vs wt7'] = sum(abs(hassm.wt_1 - hassm.wt_7))
    branchingjson['(wt1 & wt7) vs wt18'] = sum(abs((hassm.wt_1 * hassm.wt_7) - hassm.wt_18))
    branchingjson['(wt1 & wt7 & wt18) vs wt19'] = sum((hassm.wt_1 * hassm.wt_7 * hassm.wt_18) != hassm.wt_19)
    branchingjson['(mut11.1 & mut11.2 & mut15) vs wt19'] = sum(hassm.wt_19 != (hassm.mut_11_1 * hassm.mut_11_2 * hassm.mut_15))
    branchingjson['(mut11.1 & mut11.2 & mut15) vs (wt1 & wt7 & wt18)'] = sum( (hassm.mut_11_1 * hassm.mut_11_2 * hassm.mut_15) != (hassm.wt_1 * hassm.wt_7 * hassm.wt_18))
    branchingjson['All'] = sum(hassm.wt_1 * hassm.wt_7 * hassm.wt_18 * hassm.wt_19 * hassm.mut_11_1 * hassm.mut_11_2 * hassm.mut_15)
    with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/axillary_plots/branching_sm_stats.json', 'w') as fout:
        json.dump(branchingjson, fout, indent=2)

    # newick_string = "(((wt_1:24,wt_7:24):6,wt_18:6):9,wt_19:48,((mut_11_1:17.5,mut_11_2:17.5):18.5,mut_15:18.5):15)"
    # tree = Phylo.read(StringIO(newick_string), "newick")
    # tree.ladderize()
    # fig = plt.figure(figsize=[3,3], dpi=300)
    # ax = fig.add_subplot()
    # Phylo.draw(tree, axes=ax)
    # ax = cleanax(ax)
    # ax.set_xlabel("Genomic distance")
    # ax.set_ylabel("")
    # ax.spines['left'].set_visible(False)
    # plt.tick_params(axis='y', left=False, labelleft=False)

    # </editor-fold>


    # <editor-fold desc="Grouping of SMs based on branching">
    # Generate the Figure 5 from Tomimoto and Satake, 2023
    l1cnts = defaultdict(int)
    l2cnts = defaultdict(int)
    leafcnts = defaultdict(int)
    anycnts = defaultdict(int)
    for grp in allsmmat.groupby('chromosome position alt_allele tissue'.split()):
        selb = tuple([1 if branch in grp[1].branch.to_list() else 0 for branch in branchorder])
        if selb == (0, 0, 0, 0, 1, 0, 1):
            print(grp[0], selb)
        if selb == (1, 1, 1, 1, 0, 1, 1):
            print(grp[0], selb)
        if selb == (0, 1, 1, 0, 0, 0, 0):
            print(grp[0], selb)
        if grp[0][3] == 'L1':
            l1cnts[selb] += 1
        elif grp[0][3] == 'L2':
            l2cnts[selb] += 1
        elif grp[0][3] == 'leaf':
            leafcnts[selb] += 1
    for grp in allsmmat.groupby('chromosome position alt_allele'.split()):
        selb = tuple([1 if branch in grp[1].branch.to_list() else 0 for branch in branchorder])
        if selb == (1, 0, 0, 1, 0, 0, 0):
            print(grp[0], selb)
        anycnts[selb] += 1

    bkeys = set(list(l1cnts) + list(l2cnts) + list(leafcnts) + list(anycnts))
    bkeys = sorted(sorted(bkeys), key=lambda x: sum(x))

    cnts = {'L1': l1cnts, 'L2': l2cnts, 'leaf': leafcnts, 'any': anycnts}
    fig = plt.figure(figsize=[4, 6], dpi=300)
    for i, l in enumerate('L1 L2 leaf any'.split()):
        data = cnts[l]
        ax = fig.add_subplot(4, 1, i+1)
        ax.bar(x=[str(b) for b in bkeys], height=[data[k] for k in bkeys], color=colour[l], width=0.5)
        ax.tick_params(axis='x',
                       which='both',
                       bottom=False,
                       top=False,
                       labelbottom=False)
        ax.set_ylim([0, 45])
        plt.ylabel(f'SM count in {l}')
        ax.vlines(x=6.5, ymin=0, ymax=45, color='black', linestyle='dotted')
        ax=cleanax(ax)
        ax.grid(which='major', axis='x')
        ax.set_axisbelow(True)
    ax.tick_params(axis='x',
                   which='both',
                   bottom=True,
                   top=False,
                   labelbottom=True, rotation=90)
    ax.set_xlabel('Selected branch indices')
    ax = cleanax(ax)
    plt.savefig(f'{cwd}/sm_grouping_over_branches.png')
    plt.close()





    # </editor-fold>


    # <editor-fold desc="Check consistency to the dPCR results">
    dpcrpos = (('CUR5G',8809170), ('CUR1G',11957090), ('CUR1G',38054316), ('CUR1G',40001832), ('CUR2G',15761922), ('CUR5G',8809170), ('CUR6G',5258652), ('CUR6G',12539473), ('CUR6G',15553688), ('CUR6G',15729031), ('CUR6G',19101663), ('CUR6G',20953708), ('CUR7G',3760829), ('CUR7G',10743342), ('CUR7G',14411927), ('CUR8G',2951392), ('CUR1G',22927063), ('CUR7G',722256), ('CUR2G',15095953), ('cur_utg000166l_pilon',133594))
    for pos in dpcrpos:
        print(pos, allsmmat.loc[(allsmmat.chromosome == pos[0]) & (allsmmat.position == pos[1])].shape)
        print(pos)

    # </editor-fold>


    # <editor-fold desc="Not Used: Get mutations that are different between Leaf and paired L2">
    branches = ['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']
    for b in branches:
        pos = defaultdict(set)
        for tissue in 'leaf L2'.split():
            try:
                pos[tissue] = set(allsmmat.loc[allsmmat[f'{b}_{tissue}'] == 1].index.values)
            except KeyError:
                pass
        # print(b, 'leaf', pos['leaf'].difference(pos['L2']))
        # print(b, 'L2', pos['L2'].difference(pos['leaf']))
        onlyleaf = pos['leaf'].difference(pos['L2'])
        for s in onlyleaf:
            if sum(allsmmat.loc[s, allsmmat.columns != f'{b}_leaf']) > 0:
                print(b, 'leaf', s)
        onlyfruit = pos['L2'].difference(pos['leaf'])
        for s in onlyfruit:
            if sum(allsmmat.loc[s, allsmmat.columns != f'{b}_L2']) > 0:
                print(b, 'leaf', s)
    ## This found two positions
    ## mut_11_1 leaf ('CUR2G', 367749, '+TA'): Present in all samples expect mut_11_1_l2
    ## mut_15 leaf ('CUR3G', 8510980, '-GAG'): Found because it is present in both L1 and L2
    # </editor-fold>

    return
# END


def expression_plots():
    """
        Collection of plots for visualising expression differences.
    """
    # <editor-fold desc="Define imports">
    import os
    from glob import glob
    from itertools import product
    import pandas as pd
    from collections import deque, Counter
    import pybedtools as bt
    import seaborn as sns
    from matplotlib import pyplot as plt
    from hometools.strings import grepl
    from hometools.hometools import revcomp
    import json
    import pickle

    # </editor-fold>


    # <editor-fold desc="Define defaults and constants">
    basedict = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
    sdict = dict(zip(('WT_1', 'wt7', 'wt18', 'WT_19', 'MUT_11_1', 'mut11_2', 'mut4', 'MUT_15'), ('wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_4', 'mut_15')))
    rnadir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/scrna_clusters/'
    isodir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/'
    scdir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/sahu_analysis/analysis/'
    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'

    SMALL_SIZE = 7
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 8
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=BIGGER_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    plt.rc('font', family='Arial')  # fontsize of the figure title
    # </editor-fold>

    muts = pd.read_table(f'{cwd}/all_sm_in_all_samples.manually_selected.cleaned.csv', sep=',')


    # <editor-fold desc="Plot the read counts distribution of scRNA and scISO-seq reads">
    sdict2 = {'WT_1': 'wt1', 'WT_19': 'wt19', 'MUT_11_1': 'mut_11_1', 'MUT_15': 'mut_15'}
    for i, (k, v) in enumerate(sdict2.items()):
        bcrc = pd.read_table(f'{scdir}/{v}_bc_readcnts.txt', header=None, delimiter=' ')
        df = pd.read_table(f'{isodir}/{k}/cnt.txt', header=None, delim_whitespace=True)
        df[1] = df[1].str.replace('CB:Z:', '')
        df[1] = list(map(revcomp, df[1]))           # CB in Iso-Seq is reverse complemented
        df[0] = df[0].astype(int)
        isoillrc = pd.merge(bcrc, df, how='inner', on=[1])
        isoillrc.columns = "scRNA-rc BC scIso-rc".split()
        isoillrc['branch'] = k
        g = sns.jointplot(data=isoillrc, x="scRNA-rc", y="scIso-rc", kind="hex", height=3)
        g.set_axis_labels("scRNA read count", "scIso-Seq read count")
        g.fig.suptitle(f"Read count in {bnamedict[k.lower()]}")
        plt.tight_layout()
        g.savefig(f"{cwd}/{k}_bc_read_dist.png", dpi=300, transparent=True)
        g.savefig(f"{cwd}/{k}_bc_read_dist.pdf", dpi=300, transparent=True)
        plt.close()
    # </editor-fold>


    # <editor-fold desc="Plot Alt-readcount frequency">
    df = muts.sort_values(['branch', 'chromosome', 'position'])
    # # Get gene coordinates (NOT FILTERING BASED ON GENES ANYMORE)
    # gff = bt.BedTool('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.pasa_out.sort.protein_coding.3utr.gff3').saveas()
    # genes = gff.filter(lambda x: x[2] == 'gene').saveas()
    # pos = df.copy()
    # pos['start'] = pos.position - 1001
    # pos['end'] = pos.position + 1000
    # pos = pos[['chromosome', 'start', 'end', 'position']].copy()
    # pos.drop_duplicates(inplace=True)
    # posbt = bt.BedTool.from_dataframe(pos)
    # genepos = posbt.intersect(genes, u=True).to_dataframe()
    # genepos.columns = ['chromosome', 'start', 'end', 'position']
    # df = df.merge(genepos, on='chromosome position'.split(), how='inner')

    df = df.loc[[True if alt[0] not in '+-' else False for alt in df.alt_allele]]
    df = df.loc[df.branch.isin(['wt_1', 'wt_19', 'mut_11_1', 'mut_15'])]
    pos = set(zip(df.chromosome, df.position, df.alt_allele))
    clstids = list(map(lambda x: f'{x[0]}_{x[1]}', product(['wt_1', 'wt_19', 'mut_11_1', 'mut_15'], [f'clstrs_{i}' for i in range(0, 15)])))

    # <editor-fold desc="Read AF in rnaseq reads">
    rnarc = pd.DataFrame()
    for bname in ('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15'):
        clsrc = sorted(glob(f'{rnadir}/{bname}/clstrs_*.rna.rc.txt'))
        for cls in clsrc:
            clsdf = defaultdict(dict)
            with open(cls, 'r') as fin:
                for line in fin:
                    line = line.strip().split()
                    if len(line) == 0:
                        continue
                    alts = deque()
                    for p in pos:
                        if (line[0], int(line[1])) == (p[0], p[1]):
                            alts.append(p[2])
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
            clsname = cls.split('/')[-1].split(".")[0]
            # clsname = '_'.join(['clstrs', str(int(clsname.split('_')[1]) - 1)])
            clsname = '_'.join(['clstrs', str(int(clsname.split('_')[1]))])
            clsdf['cls'] = clsname
            clsdf['branch'] = sdict[bname]
            rnarc = pd.concat([rnarc, clsdf])
    rnarc.reset_index(level=[0, 1, 2], inplace=True)
    rnarc.columns = ['chromosome', 'position', 'alt_allele'] + list(rnarc.columns[3:])
    rnatrans = rnarc.loc[rnarc.rc > 0]
    rnarc = rnarc.loc[~(rnarc.af == 0)]
    rnarc['p'] = rnarc.chromosome.astype(str) + '_' + rnarc.position.astype(str) + '_' + rnarc.alt_allele.astype(str)
    rnarc['c'] = rnarc.branch.astype(str) + '_' + rnarc.cls.astype(str)
    # clstids = list(map(lambda x: f'{x[0]}_{x[1]}', product(['wt_1', 'wt_19', 'mut_11_1', 'mut_15'], [f'clstrs_{i}' for i in range(1, 13)])))
    rnapos = pd.DataFrame(product(rnarc.p.unique(), clstids))
    rnapos.columns = ['p', 'c']
    rnapos = rnapos.merge(rnarc, how='outer', on=['p', 'c'])
    rnapos.loc[rnapos.rc.isna(), ['rc', 'ac', 'af']] = 0
    # </editor-fold>


    # <editor-fold desc="Read AF in isoseq reads">
    isorc = pd.DataFrame()
    for bname in ('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15'):
        clsrc = sorted(glob(f'{isodir}/{bname}/clstrs_*.iso_seq.rc.txt'))
        for cls in clsrc:
            clsdf = defaultdict(dict)
            # with open('{}/{}/{}'.format(cwd, bname, cls), 'r') as fin:
            with open(cls, 'r') as fin:
                for line in fin:
                    line = line.strip().split()
                    if len(line) == 0: continue
                    alts = deque()
                    for p in pos:
                        if (line[0], int(line[1])) == (p[0], p[1]):
                            alts.append(p[2])
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
            clsname = cls.split('/')[-1].split(".")[0]
            # clsname = '_'.join(['clstrs', str(int(clsname.split('_')[1]) - 1)])
            clsname = '_'.join(['clstrs', str(int(clsname.split('_')[1]))])
            clsdf['cls'] = clsname
            clsdf['branch'] = sdict[bname]
            isorc = pd.concat([isorc, clsdf])
    isorc.reset_index(level=[0, 1, 2], inplace=True)
    isorc.columns = ['chromosome', 'position', 'alt_allele'] + list(isorc.columns[3:])
    isotrans = isorc.loc[isorc.rc > 0]
    isorc = isorc.loc[~(isorc.af == 0)]
    # isorc = isorc.loc[~(isorc.ac < 2)]
    isorc['p'] = isorc.chromosome.astype(str) + '_' + isorc.position.astype(str) + '_' + isorc.alt_allele.astype(str)
    isorc['c'] = isorc.branch.astype(str) + '_' + isorc.cls.astype(str)
    # clstids = list(map(lambda x: f'{x[0]}_{x[1]}', product(['wt_1', 'wt_19', 'mut_11_1', 'mut_15'], [f'clstrs_{i}' for i in range(1, 13)])))
    isopos = pd.DataFrame(product(isorc.p.unique(), clstids))
    isopos.columns = ['p', 'c']
    isopos = isopos.merge(isorc, how='left', on=['p', 'c'])
    isopos.loc[isopos.rc.isna(), ['rc', 'ac', 'af']] = 0
    # </editor-fold>


    # <editor-fold desc="Filter and cluster somatic mutations">
    # Remove position genotyped because of mismapping reads
    rnapos = rnapos.loc[rnapos.p != 'CUR7G_4335333_A']
    isopos = isopos.loc[isopos.p != 'CUR7G_4335333_A']
    mutsfilt = df.copy()
    mutsfilt['p'] = mutsfilt.chromosome.astype(str) + '_' + mutsfilt.position.astype(str) + '_' + mutsfilt.alt_allele.astype(str)
    mutsfilt = mutsfilt.loc[mutsfilt.p.isin(set(isopos.p) | set(rnapos.p))]
    mutsfilt.index = mutsfilt.p
    mutsfilt.drop('chromosome position alt_allele p'.split(), inplace=True, axis=1)
    mutsfilt = mutsfilt.loc[mutsfilt.branch.isin(['wt_1', 'wt_19', 'mut_11_1', 'mut_15'])]
    mutsfilt['sample'] = mutsfilt.branch.astype(str) + '_' + mutsfilt.tissue.astype(str)
    mutsfilt['value'] = 1
    mutsfilt = mutsfilt.pivot(columns='sample', values='value')
    mutsfilt.fillna(0, inplace=True)
    samids = list(map(lambda x: f'{x[0]}_{x[1]}', product(['wt_1', 'wt_19', 'mut_11_1', 'mut_15'], 'leaf L2'.split()))) + 'wt_1_L1 wt_19_L1 mut_11_1_L1 mut_15_L1'.split()
    for s in samids:
        if s not in mutsfilt.columns:
            mutsfilt[s] = 0.0
    mutsfilt = mutsfilt.loc[:, samids]
    mutsfilt = mutsfilt.loc[(mutsfilt != 0).any(axis=1)]

    # Cluster SMs together
    distance_matrix = pdist(mutsfilt)
    linkage_matrix = linkage(distance_matrix, method='single')
    leaves = leaves_list(optimal_leaf_ordering(linkage_matrix, distance_matrix))
    ind = mutsfilt.index.values
    indorder = deque()
    for i in leaves:
        if len(grepl('L1', mutsfilt.columns[np.where(mutsfilt.loc[ind[i]] == 1)[0]].to_list())) > 0:
            indorder.append(i)
            continue
    for i in leaves:
        if len(grepl('L2', mutsfilt.columns[np.where(mutsfilt.loc[ind[i]] == 1)[0]].to_list())) > 0:
            indorder.append(i)
            continue
        if len(grepl('leaf', mutsfilt.columns[np.where(mutsfilt.loc[ind[i]] == 1)[0]].to_list())) > 0:
            indorder.append(i)
            continue

    # mutsfilt = mutsfilt.iloc[leaves_list(optimal_leaf_ordering(linkage_matrix, distance_matrix))]
    mutsfilt = mutsfilt.iloc[indorder]
    mutsfilt = mutsfilt.T
    # Save positions of SM SNPs which are expressed in scRNA data
    snps = [i.split('_') for i in mutsfilt.columns]
    snps = pd.DataFrame(snps)
    snps[1] = snps[1].astype(int)
    snps[3] = snps[2]
    snps[2] = snps[1]
    snps[1] -= 1
    snps.sort_values([0, 1], inplace=True)
    snps.to_csv(f'{cwd}/all_sm_snps_expressed.txt', index=False, header=False, sep='\t')
    # </editor-fold>


    # <editor-fold desc="Fig1 without splitting samples based on clusters">
    cmap = sns.blend_palette(["white", colour['theme1']], as_cmap=True)
    def getmutlist(df):
        hassm = defaultdict(dict)
        for grp in df.groupby('p branch'.split()):
            if grp[1].shape[0] > 0:
                hassm[grp[0][0]][grp[0][1]] = 1
            else:
                hassm[grp[0][0]][grp[0][1]] = 1
        rnahm = pd.DataFrame.from_dict(hassm)
        tofillindex = [b for b in 'wt_1 wt_19 mut_11_1 mut_15'.split() if b not in rnahm.index.values]
        rnahm = pd.concat([rnahm, pd.DataFrame([], index=tofillindex)])
        rnahm.fillna(0, inplace=True)
        rnahm = rnahm.loc['wt_1 wt_19 mut_11_1 mut_15'.split()]
        for c in mutsfilt.columns:
            if c not in rnahm.columns:
                rnahm[c] = 0
        rnahm = rnahm.loc[:, mutsfilt.columns]
        return rnahm
    # END

    fig = plt.figure(figsize=[4, 4], dpi=300)
    ax = plt.subplot2grid((4, 1), (0, 0), rowspan=2, colspan=1, fig=fig)
    ax = sns.heatmap(mutsfilt, linewidths=0.1, linecolor='lightgrey', cmap=cmap, xticklabels=False, yticklabels=mutsfilt.index, cbar_kws={'label': 'Somatic mutation', 'fraction': 0.05}, cbar=False, ax=ax)
    ax.hlines([8], *ax.get_xlim(), color='k')
    ax.set_ylabel('Somatic mutation')
    ax.set_xlabel('')

    rnahm = getmutlist(rnapos)
    ax = plt.subplot2grid((4, 1), (2, 0), rowspan=1, colspan=1)
    ax = sns.heatmap(rnahm, linewidths=0.1, linecolor='lightgrey', cmap=cmap, xticklabels=False, yticklabels=rnahm.index, cbar_kws={'label': 'scRNA-Seq', 'fraction': 0.05}, cbar=False, ax=ax)
    ax.set_ylabel('scRNA-Seq')
    ax.set_xlabel('')

    isohm = getmutlist(isopos)
    ax = plt.subplot2grid((4, 1), (3, 0), rowspan=1, colspan=1)
    ax = sns.heatmap(isohm, linewidths=0.1, linecolor='lightgrey', cmap=cmap, yticklabels=isohm.index, cbar_kws={'label': 'scIso-Seq', 'fraction': 0.05}, cbar=False, ax=ax)
    ax.set_ylabel('scIso-Seq')
    ax.set_xlabel('')
    plt.tight_layout()
    plt.savefig('/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_af_in_iso_seq_no_clusters.png')
    plt.close()

    with open(f'{cwd}/rna.alt_readcount.pickle', 'wb') as f:
        pickle.dump({'mutsfilt': mutsfilt, 'rnahm': rnahm, 'isohm': isohm, 'rnatrans': rnatrans, 'isotrans': isotrans, 'mutsdf': df}, f)
    # </editor-fold>


    # <editor-fold desc="Secondary figure splitting samples based on clusters">
    fig = plt.figure(figsize=[5, 12], dpi=300)
    ax = plt.subplot2grid((7, 1), (0, 0), rowspan=1, colspan=1, fig=fig)
    ax = sns.heatmap(mutsfilt, linewidths=0.1, linecolor='lightgrey', cmap=cmap, xticklabels=False, yticklabels=mutsfilt.index, cbar_kws={'label': 'Somatic mutation', 'fraction': 0.05}, ax=ax)
    ax.hlines([8], *ax.get_xlim(), color='k')
    # ax.vlines([7], *ax.get_ylim(), color='k')
    ax.set_ylabel('')
    ax.set_xlabel('')

    ax = plt.subplot2grid((7, 1), (1, 0), rowspan=3, colspan=1)
    rnahm = rnapos.pivot(index='c', columns='p')['af']
    rnahm = rnahm.loc[clstids]
    for c in mutsfilt.columns:
        if c not in rnahm.columns:
            rnahm[c] = 0
    rnahm = rnahm.loc[:, mutsfilt.columns]
    ax = sns.heatmap(rnahm, linewidths=0.1, linecolor='lightgrey', cmap=cmap, xticklabels=False, yticklabels=rnahm.index, cbar_kws={'label': 'scRNA-Seq AF', 'fraction': 0.05}, ax=ax)
    ax.hlines([14, 28, 42], *ax.get_xlim(), color='k')
    # ax.vlines([7], *ax.get_ylim(), color='k')
    ax.set_ylabel('')
    ax.set_xlabel('')

    ax = plt.subplot2grid((7, 1), (4, 0), rowspan=3, colspan=1)
    isohm = isopos.pivot(index='c', columns='p')['af']
    isohm = isohm.loc[clstids]
    for c in mutsfilt.columns:
        if c not in isohm.columns:
            isohm[c] = 0
    isohm = isohm.loc[:, mutsfilt.columns]
    ax = sns.heatmap(isohm, linewidths=0.1, linecolor='lightgrey', cmap=cmap, yticklabels=isohm.index, cbar_kws={'label': 'scIsoSeq AF', 'fraction': 0.05}, ax=ax, vmin=0)
    ax.hlines([14, 28, 42], *ax.get_xlim(), color='k')
    # ax.vlines([7], *ax.get_ylim(), color='k')
    ax.set_ylabel('')
    ax.set_xlabel('')
    plt.tight_layout()
    # plt.savefig('/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_af_in_iso_seq_clusters.pdf')
    plt.savefig('/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_af_in_iso_seq_clusters.png')
    plt.close()
    # </editor-fold>


    ## Positions plotted:
    # CUR7G 4335333 Intergene
    # CUR1G 11957090    Exon
    # CUR6G 20953708    Intron
    # CUR6G 3442340 Intron
    # CUR5G 8809170 Intergene
    # CUR6G 19101663    Intergene
    # CUR5G 14960198    Intron

    # </editor-fold>

    return
# END


def gene_conversion_spectra():
    """
    Read and generate summary plots from identified gene conversions
    """
    import pandas as pd
    import pybedtools as bt

    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
    fin = f'{cwd}candidate_gene_conversion.all_sample.selected_positions.csv'
    df = pd.read_table(fin, header=None)
    df.columns = ['chromosome', 'pos1', 'pos2', 'mismatch_rc', 'match_rc', 'change_type', 'layer', 'snp_dist',	'ref_allele', 'alt_allele', 'l1_af', 'l1_rc', 'l2_af', 'l2_rc', 'l3_af',	'l3_rc', 'afdiff', 'branch', 'selected', 'position', 'ploidy', 'remark']
    # Reformat the dataframe
    dffilter = pd.DataFrame()
    for grp in df.groupby('position'):
        if grp[0] == 'CUR6G 21020736 ': continue
        for grp2 in grp[1].groupby('branch'):
            outdict = {}
            outdict['chromosome'] = [grp[0].split()[0]]
            outdict['position'] = [int(grp[0].split()[1])]
            outdict['branch'] = [grp2[0]]
            if len(set(grp2[1]['l1_af'])) != 1:
                if outdict['position'] == [30617265]:
                    garb = grp2[1].loc[grp2[1].ref_allele == 'T']
                    for k in ['layer', 'ref_allele', 'alt_allele', 'l1_af', 'l1_rc', 'l2_af', 'l2_rc', 'l3_af',	'l3_rc', 'afdiff']:
                        outdict[k] = [garb[k].unique()[0]]
                else:
                    print(f'ERROR: \n {grp}')
                    break
            else:
                for k in ['layer', 'ref_allele', 'alt_allele', 'l1_af', 'l1_rc', 'l2_af', 'l2_rc', 'l3_af',	'l3_rc', 'afdiff']:
                    outdict[k] = [grp2[1][k].unique()[0]]
            dffilter = pd.concat([dffilter, pd.DataFrame(outdict)])
    dffilter.to_csv(f'{cwd}candidate_gene_conversion.all_sample.filtered.tsv', index=False, header=True, sep='\t')

    # SMs overlapping genes, TEs, and intergenic, annotations
    gff = bt.BedTool('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.pasa_out.sort.protein_coding.3utr.gff3').saveas()
    te_rep = bt.BedTool('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/repeat/RepeatMasker/cur.genome.v1.fasta.ann.gff3').saveas()

    genes = gff.filter(lambda x: x[2]=='gene').saveas()
    cds = gff.filter(lambda x: x[2]=='CDS').saveas()
    intron = genes.subtract(gff.filter(lambda x: x[2] not in {'gene', 'mRNA'})).saveas()
    utr = gff.filter(lambda x: x[2] in {'five_prime_UTR', 'three_prime_UTR'}).saveas()
    colnames = ['chromosome', 'start', 'end', 'layer']
    pos = dffilter.copy()
    pos['start'] = pos.position - 1
    pos['end'] = pos.position
    pos = pos[colnames].copy()
    pos.drop_duplicates(inplace=True)
    posbt = bt.BedTool.from_dataframe(pos)
    genepos = posbt.intersect(genes, u=True).to_dataframe()
    cdspos = posbt.intersect(cds, u=True).to_dataframe()
    if cdspos.shape[0] > 0:
        cdspos.columns = pos.columns
        cdspos['anno'] = 'cds'
    else:
        cdspos = pd.DataFrame(columns=pos.columns.to_list() + ['anno'])
    intronpos = posbt.intersect(intron, u=True).to_dataframe()
    if intronpos.shape[0] > 0:
        intronpos.columns = pos.columns
        intronpos['anno'] = 'intron'
    else:
        intronpos = pd.DataFrame(columns=pos.columns.to_list() + ['anno'])
    utrpos = posbt.intersect(utr, u=True).to_dataframe()
    if utrpos.shape[0] > 0:
        utrpos.columns = pos.columns
        utrpos['anno'] = 'utr'
    else:
        utrpos = pd.DataFrame(columns=pos.columns.to_list() + ['anno'])
    tepos = posbt.intersect(te_rep, u=True).to_dataframe()
    if tepos.shape[0] > 0:
        tepos.columns = pos.columns
        tepos['anno'] = 'te'

    anno = np.array(['intergenic']*pos.shape[0])
    g = pos.merge(tepos, how='left', on=colnames)
    anno[g.anno == 'te'] = 'te'
    g = pos.merge(intronpos, how='left', on=colnames)
    anno[g.anno == 'intron'] = 'intron'
    g = pos.merge(utrpos, how='left', on=colnames)
    anno[g.anno == 'utr'] = 'utr'
    g = pos.merge(cdspos, how='left', on=colnames)
    anno[g.anno == 'cds'] = 'cds'
    pos['anno'] = anno

    # Out of the six gene conversions, four overlapped TEs and two in intergenic
    return
# END


def axillary_stats_plots():
    """
    Functions and commands to get extra plots for manuscript
    """

    # <editor-fold desc="Define Import">
    from gzip import open as gzopen
    from collections import deque
    from matplotlib import pyplot as plt
    import json
    from collections import defaultdict
    from glob2 import glob
    import pandas as pd
    from hometools.hometools import undict

    # </editor-fold>

    # <editor-fold desc="Define Defaults">
    outdir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/axillary_plots/'
    # </editor-fold>

    def get_fastq_readsizes(f):
        sizes = deque()
        with gzopen(f, 'r') as fin:
            for i, line in enumerate(fin):
                if i % 4 == 1:
                    sizes.append(len(line.strip()))
        return sizes
    # END

    # <editor-fold desc="Get read count, read length distribution, and total number of bases for HiFi reads">
    hififin = '/srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4784/4784_A_run521_HIFI.fastq.gz'
    sizes = get_fastq_readsizes(hififin)
    hifijsondata = dict()
    hifijsondata['Number of reads'] = len(sizes)
    hifijsondata['Total read length'] = sum(sizes)
    plt.hist(sizes, bins=100)
    plt.ylabel("Number of reads")
    plt.xlabel("Read length")
    plt.tight_layout()
    plt.savefig(f'{outdir}/hifi_read_lenght.png', dpi=300)
    plt.close()
    with open(f'{outdir}/hifi_reads.json', 'w') as fout:
        json.dump(hifijsondata, fout, indent=2)
    # </editor-fold>


    # <editor-fold desc="Get the sequencing read stats for the layer specific sample sequencing">
    indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/fruit_layer_specific/'
    branches = 'mut_11_1/  mut_11_2/  mut_15/  wt_1/  wt_18/  wt_19/  wt_7/'.replace('/', '').split()
    stats = defaultdict(dict)
    for branch in branches:
        for layer in 'l1 l2 l3'.split():
            print(branch, layer)
            stats[branch][layer] = dict()
            for i in [1, 2]:
                total_size = 0
                total_cnt = 0
                # Get read sizes for mate 1
                for f in glob(f'{indir}/{branch}/{branch}_{layer}/*{i}.fq.gz'):
                    sizes = get_fastq_readsizes(f)
                    total_size += sum(sizes)
                    total_cnt += len(sizes)
                stats[branch][layer][i] = {'size': total_size, 'count': total_cnt}
            stats[branch][layer] = {'size': stats[branch][layer][1]['size'] + stats[branch][layer][2]['size'],
                                    'count': stats[branch][layer][1]['count'] + stats[branch][layer][2]['count']}
    stats_flat = undict(stats)[0]
    df = pd.DataFrame(stats_flat)
    df = df.pivot(index=[0, 1], columns=[2]).reset_index()
    df.columns = 'branch layer sequenced_reads sequenced_bases'.split()
    df['coverage'] = df['sequenced_bases']/242500000
    df.to_csv(f'{outdir}/layer_sequencing_stats.tsv', sep='\t', index=False, header=True)
    # </editor-fold>


    # <editor-fold desc="Get the sequencing read stats for the leaf sample sequencing">
    indir1 = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/bigdata/'
    indir2 = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/rp_branch_specific/'
    outdir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/axillary_plots/'
    stats = defaultdict(dict)
    for branch in 'WT_1 WT_19 MUT_11_1 MUT_15'.split():
        for i in [1, 2]:
            total_size, total_cnt = 0, 0
            for f in glob(f'{indir1}/{branch}/{branch}*R{i}_001.fastq.gz'):
                sizes = get_fastq_readsizes(f)
                total_size += sum(sizes)
                total_cnt += len(sizes)
            stats[branch.lower()][i] = {'size': total_size, 'count': total_cnt}
    for branch in 'wt_7 wt_18 mut_4 mut_11_2'.split():
        for i in [1, 2]:
            total_size, total_cnt = 0, 0
            for f in glob(f'{indir2}/{branch.replace("_", "", 1)}/*{i}.fq.gz'):
                sizes = get_fastq_readsizes(f)
                total_size += sum(sizes)
                total_cnt += len(sizes)
            stats[branch][i] = {'size': total_size, 'count': total_cnt}

    for branch in 'mut_11_1  mut_11_2  mut_15 mut_4  wt_1  wt_18  wt_19  wt_7'.split():
        stats[branch] = {'size': stats[branch][1]['size'] + stats[branch][2]['size'],
                         'count': stats[branch][1]['count'] + stats[branch][2]['count']}

    stats_flat = undict(stats)[0]
    df = pd.DataFrame(stats_flat)
    df = df.pivot(index=[0], columns=[1]).reset_index()
    df.columns = 'branch sequenced_reads sequenced_bases'.split()
    df['coverage'] = df['sequenced_bases']/242500000
    df.to_csv(f'{outdir}/leaf_sequencing_stats.tsv', sep='\t', index=False, header=True)

    # </editor-fold>


    # <editor-fold desc="Get number of reads and size for raw scRNA-seq data">
    indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scrna/bigdata/'
    outdir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/axillary_plots/'
    stats = defaultdict(dict)
    for branch in 'WT_1 WT_19 MUT_11_1 MUT_15'.split():
        for i in [1, 2]:
            total_size, total_cnt = 0, 0
            for f in glob(f'{indir}/{branch}/*R{i}_001.fastq.gz'):
                print(branch, f)
                sizes = get_fastq_readsizes(f)
                total_size += sum(sizes)
                total_cnt += len(sizes)
            stats[branch.lower()][i] = {'size': total_size, 'count': total_cnt}

    for branch in 'mut_11_1  mut_15   wt_1    wt_19'.split():
        stats[branch] = {'size': stats[branch][1]['size'] + stats[branch][2]['size'],
                         'count': stats[branch][1]['count'] + stats[branch][2]['count']}
    stats_flat = undict(stats)[0]
    df = pd.DataFrame(stats_flat)
    df = df.pivot(index=[0], columns=[1]).reset_index()
    df.columns = 'branch sequenced_reads sequenced_bases'.split()
    df['coverage'] = df['sequenced_bases']/242500000
    df.to_csv(f'{outdir}/scrna_sequencing_stats.tsv', sep='\t', index=False, header=True)
    # </editor-fold>

    # <editor-fold desc="Get number of reads and size for raw scISO-seq data">
        # done using bash commands at SH/iso_seq_analysis.sh:9
    # </editor-fold>

    return
# END


def get_overlap_of_sm_and_flowering_gene():
    """
    Have a list of Currot genes that are homologous to A.thaliana flowering time genes (
        /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/pheno_gene/cur_athal_ortho_gene_annotations.txt
    ).
    Use this to check whether the somatic mutations are on/near these flowering
    time associated genes

    Initial test (with the SMs in leaves and only MUT_11_1 SMs) was done in get_candidate_for_mutant_phenotype.py
    """
    import pandas as pd
    import pyranges as pr

    # Reads list of flower genes and list of somatic mutations and create pyranges object
    flodf = pd.read_table("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/pheno_gene/cur_athal_ortho_gene_annotations.txt")
    flodf.rename(columns={'chrom': 'Chromosome', 'start': 'Start', 'end': 'End'}, inplace=True)
    flopr = pr.PyRanges(flodf)
    smdf = pd.read_table("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.filtered.txt")
    # Set start and end as 10kb up/down-stream from the mutation loci
    smdf['Start'] = smdf['position'] - 10000
    smdf.loc[smdf.Start < 0, 'Start'] = 1
    smdf['End'] = smdf['position'] + 10000
    smdf.rename(columns={'chromosome': 'Chromosome'}, inplace=True)
    smpr = pr.PyRanges(smdf)

    # Get overlap between the SM and flogen pyrange object
    smflo = smpr.overlap(flopr)
    flosm = flopr.overlap(smpr)
    # print(smflo)
    # +--------------+------------+------------+--------------+-------+
    # | Chromosome   |   position | branch     | ref_allele   | +7    |
    # | (category)   |    (int64) | (object)   | (object)     | ...   |
    # |--------------+------------+------------+--------------+-------|
    # | CUR1G        |   17488551 | mut_11_1   | A            | ...   |
    # | CUR6G        |    3442340 | wt_1       | A            | ...   |
    # | CUR6G        |    3447331 | mut_11_1   | T            | ...   |
    # | CUR8G        |    3674930 | mut_11_1   | A            | ...   |
    # +--------------+------------+------------+--------------+-------+
    # print(flosm)
    # +---------------+--------------+-----------+-----------+-------+
    # | Currot_gene   | Chromosome   |     Start |       End | +6    |
    # | (object)      | (category)   |   (int32) |   (int32) | ...   |
    # |---------------+--------------+-----------+-----------+-------|
    # | Gene.7488     | CUR1G        |  17482277 |  17501672 | ...   |
    # | Gene.35083    | CUR6G        |   3441508 |   3443149 | ...   |
    # | Gene.35083    | CUR6G        |   3441508 |   3443149 | ...   |
    # | Gene.46239    | CUR8G        |   3668351 |   3680207 | ...   |
    # +---------------+--------------+-----------+-----------+-------+

    # Get SMs present in all MUT or all WT
    mtpos = pd.DataFrame()
    wtpos = pd.DataFrame()
    for g in smdf.groupby(['Chromosome', 'position']):
        if set(g[1].branch.unique()) == {'mut_11_1', 'mut_11_2', 'mut_15'}:
            mtpos = pd.concat([mtpos, g[1]])
            continue
        if set(g[1].branch.unique()) == {'wt_1', 'wt_18', 'wt_7'}:  # Using wt_19 resulted in no position being selected
            wtpos = pd.concat([wtpos, g[1]])
            continue
    ## Increase the window for selection of genes
    mtpos['Start'] = mtpos['position'] - 50000
    mtpos.loc[mtpos.Start < 0, 'Start'] = 1
    mtpos['End'] = mtpos['position'] + 50000
    wtpos['Start'] = wtpos['position'] - 50000
    wtpos.loc[wtpos.Start < 0, 'Start'] = 1
    wtpos['End'] = wtpos['position'] + 50000

    mtpr = pr.PyRanges(mtpos)
    wtpr = pr.PyRanges(wtpos)
    # Get overlap between the mt/wt SM and flogen pyrange object
    mtflo = mtpr.overlap(flopr)
    wtflo = wtpr.overlap(flopr)
    flomt = flopr.overlap(mtpr)
    flowt = flopr.overlap(wtpr)
    # mtflo
    # Out[28]: Empty PyRanges
    # wtflo
    # Out[29]: Empty PyRanges
    # flomt
    # Out[30]: Empty PyRanges
    # flowt
    # Out[31]: Empty PyRanges


    return
# END


def miscelleneous():
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

    return
# END


if __name__ == '__main__':
    """
    Generate coverage histogram from samtools depth output file  
    """
    parser = argparse.ArgumentParser("Generate coverage histogram from samtools depth output file")
    parser.add_argument('f', help='path to samtools depth output', type=argparse.FileType('r'))
    parser.add_argument('-m', help='minimum read coverage cutoff', default=0, type=int)
    parser.add_argument('-M', help='maximum read coverage cutoff', default=1000, type=int)
    parser.add_argument('-s', help='sample IDs separated by \':\'', type=str)
    parser.add_argument('-p', dest='p', help='prefix for output file', type=str, default='')

    args = parser.parse_args()

    # with open(args.f.name, 'r') as f:
    count = 0
    with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/read_depth_q30_Q10.txt', 'r') as f:
        line = f.readline().strip().split()
        freqdict = {}
        for i in range(len(line) - 2):
            freqdict[i] = {str(j):0 for j in range(args.m, args.M+1)}
            try:
                freqdict[i][line[2+i]] += 1
            except KeyError as e:
                pass

        for line in f:
            count += 1
            line = line.strip().split()
            for i in range(len(line) - 2):
                try:
                    freqdict[i][line[2+i]] += 1
                except KeyError as e:
                    pass
            if count == 10000000:
                break
            if count % 1000000 == 0:
                print(count)

    if args.s is not None:
        samples = args.s.split(':')
    else:
        samples = ['sample_' + str(i) for i in range(len(line) - 2)]

    from matplotlib import pyplot as plt

    fig = plt.figure()
    for i in range(len(samples)):
        ax = fig.add_subplot(4, 1, i+1)
        ax.bar(x= list(range(args.m, args.M+1)),
                height=[freqdict[i][str(j)] for j in range(args.m, args.M+1)])
        ax.set_title(samples[i])
    plt.savefig(args.p + 'read_coverage.pdf')
