import argparse


# Set colours
colour = {('L1', 'SNP'): '#1f77b4',
          ('L1', 'Indel'): '#84b4d6',
          ('L2', 'SNP'): '#ff7f0e',
          ('L2', 'Indel'): '#ffb97a',
          ('shared', 'SNP'): '#9467bd',
          ('shared', 'Indel'): '#c4abdb',
          ('fruit', 'SNP'): '#003399',
          ('fruit', 'Indel'): '#6699ff',
          ('leaf', 'SNP'): '#2ca02c',
          ('leaf', 'Indel'): '#8bcb8b',
          'L1': '#1f77b4',
          'L2': '#ff7f0e',
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


def mutation_spectra():
    """
    Plots describing mutations spectra of somatic mutations identified from the
    layer_specific sequencing
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
    # </editor-fold>


    # <editor-fold desc="Define default values">
    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
    refcur = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
    syriout = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out'
    branches = ['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']
    branchorder = ['wt_7', 'wt_1', 'wt_18', 'wt_19', 'mut_15', 'mut_11_1', 'mut_11_2']
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
    # filtering positions that are present in all but one sample as they are probably gene conversions
    sharedsm = sharedsm.loc[sharedsm.position != 367749]
    sharedsm = sharedsm.loc[sharedsm.position != 16736100]

    data = pd.concat([layersp, sharedsm])
    data['type'] = 'SNP'
    data.loc[[a[0] in '+-' for a in data.alt_allele], 'type'] = 'Indel'
    datafilter = data.loc[:, ['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele', 'read_count', 'allele_freq', 'Layer', 'type']].drop_duplicates(subset=['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele', 'Layer', 'type']).copy()    # Filtered dataframe to be used for generating statistics
    datafilter.sort_values('chromosome position'.split(), inplace=True)
    datafilter.to_csv('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.filtered.txt', index=False, header=True, sep='\t')

    # </editor-fold>


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
    layerjsondata["Number of SMs per branch"] = {k: v for k, v in datafilter.branch.value_counts().items()}
    # Is the statistically significance? No
    g = datafilter.branch.value_counts()
    z = (g - np.mean(g))/np.std(g)  # Z-transform shows that No value is outside the 95% range
    layerjsondata["Dimeric indels"] = {i: Counter([grp[0][2] for grp in datafilter.loc[datafilter.type == 'Indel'].groupby('chromosome position alt_allele'.split())])[i] for i in ['-AT', '-TA', '+AT', '+TA']}
    layerjsondata["Dimeric indels total count"] = sum([Counter([grp[0][2] for grp in datafilter.loc[datafilter.type == 'Indel'].groupby('chromosome position alt_allele'.split())])[i] for i in ['-AT', '-TA', '+AT', '+TA']])
    with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/axillary_plots/layer_sm_stats.json', 'w') as fout:
        json.dump(layerjsondata, fout, indent=2)

    # </editor-fold>


    # <editor-fold desc="Number of SMs per layer">
    fig, ax = plt.subplots(figsize=[3, 3], dpi=300)
    y_bottom = [0]*3
    for t in ('SNP', 'Indel'):
        layercnts = datafilter.loc[datafilter.type == t] \
            .drop_duplicates(subset='chromosome position alt_allele'.split()) \
            .Layer.value_counts().to_dict()
        ax.bar(['L1', 'L2', 'shared'], [layercnts[l] for l in ('L1', 'L2', 'shared')], bottom=y_bottom, color=[colour[l, t] for l in ('L1', 'L2', 'shared')])
        y_bottom = [layercnts[l] for l in ('L1', 'L2', 'shared')]
        print(y_bottom)
    ax.set_xlabel("Layer")
    ax.set_ylabel("Number of SMs")
    ax = cleanax(ax)
    # ax.legend(bbox_to_anchor=(1.01, 1), frameon=False)
    plt.tight_layout(pad=0.1)
    plt.savefig(f'{cwd}/sm_counts_in_layers.png')
    plt.close()
    # </editor-fold>


    # <editor-fold desc="Get syri plot with SMs">
    df = datafilter.drop_duplicates(subset=['chromosome', 'position', 'ref_allele', 'alt_allele'])
    for grp in df.groupby('Layer'):
        pos = grp[1].copy()
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
    for tissue in 'L1 L2 shared'.split():
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
    fig = plt.figure(figsize=[2, 2], dpi=300)
    ax = fig.add_subplot()
    ax = sns.barplot(data=smdist, x='variable', y='value', hue='tissue', hue_order='L1 L2 shared'.split(), palette=colour, ax=ax)
    ax.set_ylim([0, 1.e-8])
    ax.set_xlabel('')
    ax.set_ylabel('Ratio of SMs')
    ax = cleanax(ax)
    ax.legend(frameon=False, title='')
    plt.tight_layout(pad=0.1)
    plt.savefig(f'{cwd}/sm_overlap_with_structural_annotation.png')
    plt.close()
    # </editor-fold>


    # <editor-fold desc="Get transversion vs transitions plot">
    unisnp = datafilter.loc[datafilter.type == 'SNP'].drop_duplicates(subset=['chromosome', 'position', 'ref_allele', 'alt_allele', 'Layer'])
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
    df.columns = ['SNP type', 'Layer', 'Percentage of variants']
    fig = plt.figure(figsize=[3, 3], dpi=300)
    ax = fig.add_subplot()
    # sns.barplot(data=df, x='SNP type', y='Percentage of variants', hue='Layer', palette=colour, order=['A>G', 'T>C', 'G>A', 'C>T', 'A>C', 'T>G', 'A>T', 'T>A', 'G>T', 'C>A', 'G>C', 'C>G'], ax=ax, width=0.9)
    sns.barplot(data=df, x='SNP type', y='Percentage of variants', hue='Layer', palette=colour, order=['A,T>G,C', 'C,G>T,A', 'A,T>C,G', 'A,T>T,A', 'C,G>A,T', 'C,G>G,C'], ax=ax, width=0.75)
    ax.axvline(1.5, linestyle='dashed', color='black', linewidth=0.75)
    ax.set_ylim([0, 100])
    ax.legend(frameon=False)
    ax.tick_params(axis='x', rotation=30)
    ax = cleanax(ax)
    plt.savefig(f'{cwd}/snps_mutation_spectra.png')
    # plt.savefig(f'{cwd}/snps_mutation_spectra.pdf')
    # plot2emf(f'{cwd}/snps_mutation_spectra.pdf', f'{cwd}/snps_mutation_spectra.emf')
    plt.close()
    # </editor-fold>


    # <editor-fold desc="Get distibution of SMs SNPs in genomic triplets. Are CpGs enriched?">
    snppos = datafilter.loc[datafilter.type == 'SNP'].drop_duplicates(subset='chromosome position alt_allele'.split())
    # smvars = allsmmat.reset_index()
    # smvars['type'] = ['SNP' if i[0] not in '+-' else 'indel' for i in smvars['alt_allele']]
    # snppos = smvars.loc[smvars.type == 'SNP'].copy()   # 1 based positions

    # Count the total number of
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
    tripcntabs = defaultdict()
    tripcntnorm = defaultdict()
    for l in 'L1 L2 shared'.split():
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
    tripletjson['cgsnpcnt'] = {l: sum([tripcntabs[l][c] for c in cpgkeys]) for l in 'L1 L2 shared'.split()}
    tripletjson['othersnpcnt'] = {l: sum([v for k, v in tripcntabs[l].items() if k not in cpgkeys]) for l in 'L1 L2 shared'.split()}
    with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/axillary_plots/triplet_sm_stats.json', 'w') as fout:
        json.dump(tripletjson, fout, indent=2)
    tripkeys = sorted(kmercnts, key=lambda x: sum([tripcntnorm[l][x] for l in 'L1 L2 shared'.split()]), reverse=True)
    fig = plt.figure(figsize=[6, 3], dpi=300)
    ax = fig.add_subplot()
    ybottom = np.zeros_like(tripkeys, dtype='float')
    for l in 'L1 L2 shared'.split():
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

    # </editor-fold>


    # <editor-fold desc="Plot for number of SMs per branch">

    branchscnt = {g[0]: defaultdict(int, g[1].branch.value_counts()) for g in datafilter.groupby(['Layer', 'type'])}
    plt.rc('font', size=9)
    fig, ax = plt.subplots(figsize=[6, 3], dpi=300)
    y_bottom = [0]*7
    bar_width = 0.5
    for l in ('L1', 'L2', 'shared'):
        for t in ('SNP', 'Indel'):
            y = [branchscnt[(l, t)][b] for b in branches]
            ax.bar(np.arange(1, 8), y, label=f'{l} {t}', bottom=y_bottom, color=colour[(l, t)], width=bar_width)
            y_bottom = [y_bottom[i] + y[i] for i in range(7)]
    ax.set_xlabel("Branch")
    ax.set_ylabel("Number of SMs")
    ax.set_xticks(np.arange(1, 8))
    ax.set_xticklabels(branches)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend(bbox_to_anchor=(1.01, 1), frameon=False)
    plt.tight_layout(pad=0.1)
    plt.savefig(f'{cwd}/sm_counts_in_branches.png')
    plt.close()
    # </editor-fold>


    # <editor-fold desc="Indel size distribution">
    inddf = datafilter.loc[datafilter['type'] == 'Indel'].drop_duplicates(subset=['chromosome', 'position', 'ref_allele', 'alt_allele', 'Layer'])
    indsize = dict()
    for grp in inddf.groupby('Layer'):
        # printdf(grp[0], grp[1])
        indsize[grp[0]] = defaultdict(int, Counter([len(a)-1 if '+' in a else 1-len(a) for a in grp[1].alt_allele]))
    maxx = max(map(abs, unlist([list(v.keys()) for v in indsize.values()])))
    xran = list(range(-maxx, maxx+1))

    ybottom = np.zeros_like(xran)
    fig = plt.figure(figsize=[4, 3], dpi=300)
    ax = fig.add_subplot()
    for l in 'L1 L2 shared'.split():
        lheight = np.array([indsize[l][x] for x in xran])
        ax.bar(xran, lheight, bottom=ybottom, color=colour[l], label=l)
        ybottom += lheight
    ax.set_xlabel("Indel size")
    ax.set_ylabel("Number of indels")
    plt.tight_layout()
    ax.legend(frameon=False)
    ax = cleanax(ax)
    plt.savefig(f'{cwd}/indel_size_distribution.png')
    plt.close()
    # </editor-fold>


    # <editor-fold desc="SMs overlapping genes, TEs, and intergenic, annotations">
    ## TODO: Permutation test for the distribution of ShVs?
    ## TODO: Get SM distributions between different TE/repeat types
    gff = bt.BedTool('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.pasa_out.sort.protein_coding.3utr.gff3').saveas()
    te_rep = bt.BedTool('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/repeat/RepeatMasker/cur.genome.v1.fasta.ann.gff3').saveas()

    genes = gff.filter(lambda x: x[2] == 'gene').saveas()
    cds = gff.filter(lambda x: x[2] == 'CDS').saveas()
    intron = genes.subtract(gff.filter(lambda x: x[2] not in {'gene', 'mRNA'})).saveas()
    utr = gff.filter(lambda x: x[2] in {'five_prime_UTR', 'three_prime_UTR'}).saveas()

    pos = datafilter.copy()
    pos['start'] = pos.position - 1
    pos['end'] = pos.position
    pos = pos[['chromosome', 'start', 'end', 'Layer', 'type']].copy()
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
    g = pos.merge(tepos, how='left', on=['chromosome', 'start', 'end', 'Layer', 'type'])
    anno[g.anno == 'te'] = 'te'
    g = pos.merge(intronpos, how='left', on=['chromosome', 'start', 'end', 'Layer', 'type'])
    anno[g.anno == 'intron'] = 'intron'
    g = pos.merge(utrpos, how='left', on=['chromosome', 'start', 'end', 'Layer', 'type'])
    anno[g.anno == 'utr'] = 'utr'
    g = pos.merge(cdspos, how='left', on=['chromosome', 'start', 'end', 'Layer', 'type'])
    anno[g.anno == 'cds'] = 'cds'
    pos['anno'] = anno

    regsize = dict()
    regsize['cds'] = sum([int(b[2])-int(b[1])+1 for b in cds.sort().merge()])
    regsize['utr'] = sum([int(b[2])-int(b[1])+1 for b in utr.sort().merge()])
    regsize['intron'] = sum([int(b[2])-int(b[1])+1 for b in intron.sort().merge()])
    regsize['te'] = sum([int(b[2])-int(b[1])+1 for b in te_rep.sort().merge()])
    regsize['intergenic'] = 233151598 - sum(regsize.values())


    xticks = ['cds', 'utr', 'intron', 'te', 'intergenic']
    bar_width = 0.3
    fig = plt.figure(figsize=[4, 6], dpi=300)
    ax = fig.add_subplot(3, 1, 1)

    for i, l in enumerate('L1 L2 shared'.split()):
        y_bottom = [0]*5
        for t in 'SNP Indel'.split():
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

    ax = fig.add_subplot(3, 1, 2)
    for i, l in enumerate('L1 L2 shared'.split()):
        y_bottom = [0]*5
        for t in 'SNP Indel'.split():
            cnts = pos.loc[(pos.Layer == l) & (pos.type == t), 'anno'].value_counts()
            y = [(cnts[x]*1000000)/regsize[x] if x in cnts else 0 for x in xticks]
            ax.bar(np.arange(5) - (bar_width*(1-i)), y, label=f'{l} {t}', bottom=y_bottom, color=colour[l, t], width=bar_width)
            y_bottom = y
            ax.set_xticks(np.arange(5))
            ax.set_xticklabels(['Coding', 'UTR', 'Intron', 'TE/Repeat', 'Intergenic'])
            ax.set_ylabel('Mutations per MBp')
            ax.set_xlabel('')
            ax = cleanax(ax)

    ax = fig.add_subplot(3, 1, 3)
    for i, l in enumerate('L1 L2 shared'.split()):
        y_bottom = [0]*5
        for t in 'SNP Indel'.split():
            cnts = pos.loc[(pos.Layer == l) & (pos.type == t), 'anno'].value_counts()
            cnts = (cnts/pos.Layer.value_counts()[l])*100
            y = [cnts[x] if x in cnts else 0 for x in xticks]
            ax.bar(np.arange(5) - (bar_width*(1-i)), y, label=f'{l} {t}', bottom=y_bottom, color=colour[l, t], width=bar_width)
            y_bottom = y
            ax.set_xticks(np.arange(5))
            ax.set_xticklabels(['Coding', 'UTR', 'Intron', 'TE/Repeat', 'Intergenic'])
            ax.set_ylabel('Ratio of mutations')
            ax.set_xlabel('Annotation')
            ax = cleanax(ax)

    plt.tight_layout(pad=0)
    plt.savefig(f'{cwd}/all_layer_somatic_variants.annotation_overlap.png')
    plt.close()
    # subprocess.call(f'inkscape {cwd}/all_layer_somatic_variants.annotation_overlap.svg -M {cwd}/all_layer_somatic_variants.annotation_overlap.emf', shell=True)
    # </editor-fold>


    # <editor-fold desc="Number of branches covered by each SM">
    pos = datafilter.drop_duplicates(subset=['chromosome', 'position', 'alt_allele']).copy()
    samplecnt = {grp[0]: len(set(grp[1].branch)) for grp in datafilter.groupby(['chromosome', 'position', 'alt_allele'])}
    pos['Branch Count'] = [samplecnt[row.chromosome, row.position, row.alt_allele] for row in pos.itertuples(index=False)]
    n_branch = pos.groupby(['Layer', 'type', 'Branch Count']).size().to_dict()
    # plt.rc('font', size=15)
    fig, ax = plt.subplots(figsize=[3, 3], dpi=300)
    y_bottom = [0]*7
    bar_width = 0.5
    for l in ('L1', 'L2', 'shared'):
        for t in ('SNP', 'Indel'):
            y = [n_branch[(l, t, i)] if (l, t, i) in n_branch else 0 for i in range(1, 8)]
            ax.bar(np.arange(1, 8), y, label=f'{l}_{t}', bottom=y_bottom, color=colour[(l, t)], width=bar_width)
            y_bottom = [y_bottom[i] + y[i] for i in range(7)]
    ax.set_xlabel("Number of Branch")
    ax.set_ylabel("Number of SMs")
    ax.set_xticks(np.arange(1, 8))
    ax.legend(bbox_to_anchor=(1.01, 1), frameon=False)
    ax = cleanax(ax)
    plt.savefig(f'{cwd}/number_sm_per_branch.png')
    plt.close()
    # </editor-fold>


    # <editor-fold desc="Example plot showing mutations having higher AF in L1 in all branches">
    pos = ('CUR1G', 34111279, '+AT')
    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
    branches = ['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']
    readcnt = defaultdict(dict)
    for bname in branches:
        # for l in ('l1', 'l2', 'l3'):
        for l in ('l1', 'l2'):
            with open(f'{cwd}/{bname}/{bname}_{l}/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', 'r') as fin:
                for line in fin:
                    line = line.strip().split()
                    if line[0] == pos[0]:
                        if int(line[1]) == pos[1]:
                            try:
                                index = line.index(pos[2])
                            except ValueError:
                                readcnt[bname][l] = 0
                                break
                            readcnt[bname][l] = int(line[index+1])
                            break
                    else:
                        readcnt[bname][l] = 0
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
    plt.savefig(f"{cwd}/layer_conserved_somatic_mutation_allele_frequency.png")
    plt.close()
    # </editor-fold>


    # <editor-fold desc="Grouping of SMs based on branching">
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
    bkeys = sorted(sorted(bkeys), key=lambda x: sum(x))

    cnts = {'L1': l1cnts, 'L2': l2cnts, 'shared': shcnts}
    fig = plt.figure(figsize=[4, 6], dpi=300)
    for i, l in enumerate('L1 L2 shared'.split()):
        data = cnts[l]
        ax = fig.add_subplot(3, 1, i+1)
        ax.bar(x=[str(b) for b in bkeys], height = [data[k] for k in bkeys], color=colour[l], width=0.5)
        ax.tick_params(axis='x',
                       which='both',
                       bottom=False,
                       top=False,
                       labelbottom=False)
        ax.set_ylim([0, 25])
        plt.ylabel(f'SM count in {l}')
        ax.vlines(x=6.5, ymin=0, ymax=25, color='black', linestyle='dotted')
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


    # <editor-fold desc="Other plots that does not seem too relevant in the current description">
    # Distribution of SMs over the genome
    chrsize = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.fai', header=None)
    chrsize = chrsize.iloc[:8]
    pos = datafilter.drop_duplicates(subset=['chromosome', 'position', 'alt_allele']).copy()
    pos = pos.loc[['CUR' in c for c in pos.chromosome]]
    chroms = pd.factorize(pos.chromosome)[1]
    pos.chromosome = 7 - pd.factorize(pos.chromosome)[0]
    margin = 0.2    # Change to adjust space between L1 and L2 tracks
    pos.loc[pos.Layer == 'L1', 'chromosome'] += margin
    pos.loc[pos.Layer == 'L2', 'chromosome'] -= margin
    pos['branch_type'] = ['wt' if 'wt' in b else 'mut' for b in pos.branch]

    sns.set(rc={"figure.figsize": (4, 6)}, font_scale=1.2)
    sns.set_style(style='white')
    ax = sns.scatterplot(pos, x='position', y='chromosome', hue='type', style='type', zorder=1, s=60, alpha=0.8)
    ax.set_yticklabels([''] + list(chroms)[::-1])
    ax.set_xlim([-1000000, max(chrsize[[1]].values)[0] + 1000000])
    lab = True
    for i in range(8):
        if lab:
            for j, l in enumerate(('L1', 'L2')):
                ax.axhline(i+margin-(j*2*margin), 0, chrsize.iloc[7-i, 1]/max(chrsize[[1]].values)[0], linewidth=3, zorder=0, alpha=1, color=colour[l], label=l)
            lab = False
            continue
        for j, l in enumerate(('L1', 'L2')):
            ax.axhline(i+margin-(j*2*margin), 0, chrsize.iloc[7-i, 1]/max(chrsize[[1]].values)[0], linewidth=3, zorder=0, alpha=1, color=colour[l])
    ax.legend(frameon=False)
    ax.spines[:].set_visible(False)
    ax.tick_params(left=False)
    ax.set_ylabel(None)
    ax.set_xlabel('Position (in 10Mb)')
    ax.set_ylabel('Chromosome')
    ax.set_title("Distribution of somatic mutations\n in the genome")
    plt.tight_layout(pad=0.1)
    plt.savefig(f'{cwd}/sm_dist_over_genome.png', dpi=600)
    plt.close()


    # Get maximum spanning tree and hierarchical clustering between branches (this would show
    # whether the SMs fit the morphology of the tree)
    pos = datafilter[['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele']].copy()
    pos.index = pos.branch
    pos.drop(['branch'], axis=1, inplace=True)
    pos = pos.groupby(['chromosome', 'position', 'ref_allele', 'alt_allele']).filter(lambda x: len(x)<7)
    branchpos = {grp[0]: {tuple(r) for r in grp[1].itertuples(index=False)} for grp in pos.groupby(level=0)}
    branches = ['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']
    plt.rc('font', size=15)
    fig = plt.figure(figsize=[7, 5])
    ax = fig.add_subplot(1, 2, 1)

    ## Get distances for hierarchical clustering
    AM = deque()
    for i, s in enumerate(['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']):
        for j, e in enumerate(['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']):
            if i >= j:
                continue
            print(len(branchpos[s].intersection(branchpos[e])))
            AM.append(30 - len(branchpos[s].intersection(branchpos[e])))
    AM = list(AM)
    Z = linkage(AM, method='ward')
    dendrogram(optimal_leaf_ordering(Z, AM), ax=ax, leaf_font_size=15)
    # ax = plt.gca()
    ax.spines[:].set_visible(False)
    ax.tick_params(left=False, labelleft=False)
    ax.set_xticks(labels=[branches[i] for i in leaves_list(optimal_leaf_ordering(Z, AM))], ticks=plt.xticks()[0], rotation=90)
    # ax.set_title("Clustering of branches based on SMs")
    plt.tight_layout(pad=0.1)

    ## Get adjacency matrix for minimum spanning tree
    AM2 = 30 - squareform(AM)
    np.fill_diagonal(AM2, 0)
    g = ig.Graph.Weighted_Adjacency(AM2, mode="undirected")
    g.vs['name'] = branches
    inv_weight = [1./w for w in g.es["weight"]]
    T = g.spanning_tree(weights=inv_weight)
    print(T.is_connected())

    ax = fig.add_subplot(2, 2, 2)
    ig.plot(g, target=ax, layout='davidson_harel', vertex_label=g.vs["name"], edge_width=np.log1p(g.es['weight']), vertex_size=0.2)
    ax.set_xlabel('branch connectivity')
    ax = fig.add_subplot(2, 2, 4)
    ig.plot(T, target=ax, layout='davidson_harel', vertex_label=T.vs["name"], edge_width=np.log1p(T.es['weight']), vertex_size=0.2)
    ax.set_xlabel('Maximum spanning tree')
    plt.tight_layout(h_pad=2, w_pad=3)
    plt.savefig(f'{cwd}/sm_branch_clustering.png', dpi=600)
    plt.close()

    # </editor-fold>


    # <editor-fold desc="Get correlation of mappability and SM. Ideally, this should not be different than the background distribution">
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
    Merge the list of SMs identified in leaf, layers, and all samples together
    """

    # <editor-fold desc="Define imports">
    import pandas as pd
    from matplotlib import pyplot as plt
    import numpy as np
    import seaborn as sns
    from collections import deque, defaultdict, Counter
    import pybedtools as bt
    from matplotlib_venn import venn3
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
    branches = ['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']
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


    # <editor-fold desc="Get cleaned allsmmat file">
    layerfin = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.filtered.txt'
    leaffin = "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/high_cov_mutants_sorted.all_samples.selected.txt"
    # allfin = "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_samples_candidate_high_cov_sms3.only_selected.csv"

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

    # # Read all sample SM
    # allsam = pd.read_table(allfin, header=None)
    # allsam.columns = 'chromosome position branch tissue ref_allele alt_allele read_count allele_freq selected ploidy remark'.split()
    # allsam.tissue = allsam.tissue.str.upper()

    # allsm = pd.concat([leafsm, layersm, allsam])
    allsm = pd.concat([leafsm, layersm])
    allsm.sort_values(['chromosome', 'position'], inplace=True)

    # Change the annotation at positions incorrectly aligned manually
    allsm.loc[allsm.position == 20360183, 'type'] = 'SNP'
    allsm.loc[allsm.position == 20360183, 'alt_allele'] = 'G'

    allsm.loc[allsm.position == 19296318, 'alt_allele'] = 'A'
    allsm.loc[allsm.position == 19296318, 'ref_allele'] = 'C'
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
    # Analyse manually to select positions that are FN in other tissues
    # Read the manually curated file and update the allsmat mat
    # Removed SMs at CUR2G:367749 and CUR7G:16736100
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
    # </editor-fold>


    allsmmat = pd.read_table(f'{cwd}/all_sm_in_all_samples.manually_selected.cleaned.csv')


    # <editor-fold desc="Create summary stats">
    allsmjsondata = dict()
    allsmjsondata["Number of events"] = allsmmat.shape[0]
    allsmjsondata["Total number of unique SMs"] = allsmmat.drop_duplicates(subset='chromosome position alt_allele'.split()).shape[0]
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


    # <editor-fold desc="Number of SMs per organ">
    df = allsmmat.copy()
    df['type'] = ['SNP' if row[2][0] not in '+-' else 'Indel' for row in df.itertuples(index=False)]
    fig, ax = plt.subplots(figsize=[2, 3], dpi=300)
    y_bottom = [0]*3
    for t in ('SNP', 'Indel'):
        cnts = df.loc[df.type == t] \
            .drop_duplicates(subset='chromosome position alt_allele'.split()) \
            .organ.value_counts().to_dict()
        ax.bar(['fruit', 'leaf', 'shared'], [cnts[l] for l in ('fruit', 'leaf', 'shared')], bottom=y_bottom, color=[colour[l, t] for l in ('fruit', 'leaf', 'shared')])
        y_bottom = [cnts[l] for l in ('fruit', 'leaf', 'shared')]
        print(y_bottom)
    ax.set_xlabel("Organ")
    ax.set_ylabel("Number of SMs")
    ax = cleanax(ax)
    # ax.legend(bbox_to_anchor=(1.01, 1), frameon=False)
    plt.tight_layout(pad=0.1)
    plt.savefig(f'{cwd}/sm_counts_in_layers.png')
    plt.close()
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


    # <editor-fold desc="Venn diagram showing organ specificity of SMs">
    def tobed(df):
        df = df['chromosome position'.split()].copy()
        df.columns = 'chromosome end'.split()
        df['start'] = df.end - 1
        df = df['chromosome start end'.split()]
        return bt.BedTool.from_dataframe(df)

    leafsms = allsmmat.loc[allsmmat.tissue == 'leaf'].drop_duplicates(subset='chromosome position alt_allele'.split()).copy()
    leafsms = tobed(leafsms)
    l1sms = allsmmat.loc[allsmmat.tissue == 'L1'].drop_duplicates(subset='chromosome position alt_allele'.split())
    l1sms = tobed(l1sms)
    l2sms = allsmmat.loc[allsmmat.tissue == 'L2'].drop_duplicates(subset='chromosome position alt_allele'.split())
    l2sms = tobed(l2sms)

    fig = plt.figure(figsize=(2, 3), dpi=300)
    ax = fig.add_subplot()
    leafcnt = len(leafsms.subtract(l1sms).subtract(l2sms))
    l1cnt = len(l1sms.subtract(leafsms).subtract(l2sms))
    l2cnt = len(l2sms.subtract(leafsms).subtract(l1sms))
    leafl2cnt = len(leafsms.intersect(l2sms, u=True).subtract(l1sms))
    leafl1cnt = len(leafsms.intersect(l1sms, u=True).subtract(l2sms))
    l1l2cnt = len(l1sms.intersect(l2sms, u=True).subtract(leafsms))
    l1l2leafcnt = len(l1sms.intersect(l2sms, u=True).intersect(leafsms, u=True))
    venn3(subsets=(leafcnt, l2cnt, leafl2cnt, l1cnt, leafl1cnt, l1l2cnt, l1l2leafcnt),
          set_labels='leaf L2 L1'.split(),
          set_colors=(colour['leaf'], colour['L2'], colour['L1']), alpha=1, ax=ax)
    plt.tight_layout(rect=(0, 0, 1, 1))
    plt.savefig(f'{cwd}/organ_specificity.png')
    plt.close()
    # </editor-fold>


    # <editor-fold desc="Heatmap comparing allele frquencies of somatic variations in leaf, L1 and L2">
    ## get allsmmat regions for use with BAM-readcount
    df = allsmmat['chromosome position position'.split()].drop_duplicates()
    df.to_csv(f'{cwd}/all_sm_in_all_samples.manually_selected.cleaned.regions', index=False, header=False, sep='\t')
    ## Get readcounts from the leaf BAM files, using code in layer_specific_dna_analysis_all_samples.sh
    pos = set(list(map(tuple, allsmmat['chromosome position alt_allele'.split()].to_numpy())))
    rc = dict()
    af = dict()
    basedict = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
    sdict = dict(zip(('WT_1', 'wt7', 'wt18', 'WT_19', 'MUT_11_1', 'mut11_2', 'mut4', 'MUT_15'), ('wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_4', 'mut_15')))
    ## Read read counts at the layer SM positions in leaves
    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
    for b in ('WT_1', 'wt7', 'wt18', 'WT_19', 'MUT_11_1', 'mut11_2', 'MUT_15'):
        with open(f'{cwd}/all_sm_in_all_samples_readcounts/{b}.all_sm_in_all_samples.read_count.txt', 'r') as fin:
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
            with open(f'{cwd}/all_sm_in_all_samples_readcounts/{b}.{l}.all_sm_in_all_samples.read_count.txt', 'r') as fin:
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
    allsmorgan = allsmmat['chromosome position alt_allele organ'.split()].drop_duplicates()
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


    # <editor-fold desc="AF plots of tissue specific mutations">
    ## layer 1 mutations could not be identified in leaves, further support that
    ## accurate SM identification would require homogenous cell population)
    afvalues = pd.read_table(f'{cwd}/all_sm_in_all_samples.allele_freq.txt')
    afvalues.set_index(['Unnamed: 0', 'Unnamed: 1', 'Unnamed: 2'], inplace=True)
    ldf = allsmmat.loc[allsmmat.organ == 'leaf'].drop('tissue organ'.split(), axis=1)
    ldf.set_index('chromosome position alt_allele'.split(), inplace=True)
    leafpos = set(list(map(tuple, ldf.index.values)))
    fdf = allsmmat.loc[allsmmat.organ == 'fruit']
    fdf.set_index('chromosome position alt_allele'.split(), inplace=True)
    l1pos = set(list(map(tuple, fdf.loc[fdf.tissue == 'L1'].index.values)))
    l2pos = set(list(map(tuple, fdf.loc[fdf.tissue == 'L2'].index.values)))
    l1only = l1pos.difference(l2pos)
    l2only = l2pos.difference(l1pos)
    fig = plt.figure(figsize=[3, 4.5], dpi=300)
    for i, (k, v) in enumerate({'leaf': leafpos, 'L1': l1only, 'L2': l2only}.items()):
        ax = fig.add_subplot(3, 1, i+1)
        df = ldf if k == 'leaf' else fdf
        afcnts = dict()
        for row in afvalues.itertuples(index=True):
            if row.Index in v:
                s = df.loc[row.Index].to_numpy()[0] if k == 'leaf' else df.loc[row.Index].to_numpy()[0][0]
                afcnts[row.Index] = [s] + getvalues(row, [row._fields.index(c) for c in [s, f'{s}_l1', f'{s}_l2']])
        afcnts = pd.DataFrame(afcnts).T
        afcnts.reset_index(inplace=True)
        afcnts.columns = ['chromosome', 'position', 'alt_allele', 'sample'] + ['leaf', 'L1', 'L2']
        afcnts = afcnts.melt(id_vars=['chromosome', 'position', 'alt_allele', 'sample'])
        ax = sns.lineplot(afcnts, x='variable', y='value', units='position', color='lightgrey', alpha=0.5, linestyle='dotted',  estimator=None, ax=ax, zorder=0, size=0.1, legend=False)
        ax = sns.boxplot(data=afcnts, x='variable', y='value', ax=ax, linewidth=0.5, fliersize=0.2, palette=colour, zorder=1)
        ax.set_xlabel('')
        ax.set_ylabel('Allele Frequency')
        ax.set_ylim([0, 1])
        ax = cleanax(ax)
    plt.tight_layout(pad=0.1)
    plt.savefig(f'{cwd}/af_dist_all_sm.png')
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
    #

    # <editor-fold desc="Heatmap for presence/absence of SMs in different samples with and without clustering">
    samples = sorted(['_'.join(i) for i in product(branches, 'leaf L1 L2'.split())], key=lambda x: x.rsplit('_', 1)[1])
    allsmpivot = dict()
    for grp in allsmmat.groupby(by='chromosome  position alt_allele'.split()):
        df = grp[1].copy()
        dfsamples = set(df.branch.astype(str) + '_' + df.tissue.astype(str))
        allsmpivot[grp[0]] = {s: (1 if s in dfsamples else 0) for s in samples}
    allsmpivot = pd.DataFrame(allsmpivot).T
    isSNP = ['red' if v[2][0] in '+-' else 'blue' for v in allsmpivot.index]
    cmap = sns.blend_palette(["white", colour['theme1']], as_cmap=True)
    sns.clustermap(allsmpivot, cmap=cmap, linecolor='lightgrey', dendrogram_ratio=(0.05, 0.05), linewidths=0.01, cbar_pos=None, figsize=(8, 10), row_colors=isSNP, yticklabels=False)
    plt.ylabel('')
    plt.tight_layout()
    plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_heatmap_in_all_samples.png', dpi=300)
    plt.close()

    sns.clustermap(allsmpivot, cmap=cmap, linecolor='lightgrey', col_cluster=False, dendrogram_ratio=(0.05, 0.0), linewidths=0.01, cbar_pos=None, figsize=(8, 10), row_colors=isSNP, yticklabels=False)
    plt.ylabel('')
    plt.tight_layout()
    plt.savefig(f'{cwd}/all_sm_heatmap_in_all_samples_no_col_cluster.png', dpi=300)
    plt.close()
    # </editor-fold>


    # <editor-fold desc="Heatmap and clustering of different samples based on shared SMs">
    tissues = ['leaf', 'L1', 'L2']
    btiss = ['_'.join([b, a]) for a in tissues for b in branches]
    btissm = dict()
    for c in allsmpivot.columns:
        btissm[c] = set((allsmpivot.loc[allsmpivot[c] == 1]).index.values)
    allL1 = btissm[f'{branches[0]}_L1']
    for b in branches[1:]:
        allL1 = allL1.intersection(btissm[f'{b}_L1'])
    # Remove positions present in all L1
    for b in branches:
        btissm[f'{b}_L1'] = btissm[f'{b}_L1'].difference(allL1)

    btiss_mat = np.zeros((len(btiss), len(btiss)), dtype='int')
    for i, a in enumerate(btiss):
        for j, b in enumerate(btiss):
            btiss_mat[i][j] = len(btissm[a].intersection(btissm[b]))

    sns.clustermap(btiss_mat, cmap=cmap, xticklabels=btiss, yticklabels=btiss, linecolor='lightgrey', dendrogram_ratio=(0.05, 0.05), cbar_pos=None, linewidths=0.01, figsize=(6, 6))
    plt.tight_layout()
    plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_distribution.png', dpi=300)
    plt.close()

    for i in range(21):
        btiss_mat[i, i] = 0
    g = ig.Graph.Weighted_Adjacency(btiss_mat, mode='undirected')
    g.vs['bname'] = [i.rsplit("_", 1)[0] for i in btiss]
    g.vs['tname'] = [colour[i.rsplit("_", 1)[1]] for i in btiss]
    fig = plt.figure(figsize=[8, 8])
    ax = fig.add_subplot(1, 1, 1)
    ig.plot(g, target=ax, vertex_label=g.vs["bname"], vertex_color=g.vs["tname"], edge_width=np.log1p(g.es['weight']), layout=g.layout_auto(), margin=0, bbox=[100, 100], vertex_size=0.4)
    plt.tight_layout()
    plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_sharing.png')
    plt.close()
    # </editor-fold>


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


    # <editor-fold desc="Calculate mutation load: Consider each sample as separate lineage and divide by the diploid genome size">
    fig = plt.figure(figsize=[4, 2])
    branches = ['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']
    mrates = pd.DataFrame(allsmpivot.apply(sum, axis=0)/480000000)
    mrates.columns = ['value']
    mrates['mrate'] = 'Mutation load'
    mrates['Tissue'] = [i.rsplit('_', 1)[1] for i in mrates.index.values]
    ax = fig.add_subplot(1, 2, 1)
    ax = sns.violinplot(data=mrates, x="mrate", y="value", color=colour['theme1'], linewidth=0, inner=None, ax=ax)
    plt.setp(ax.collections, alpha=.2)
    ylim = ax.get_ylim()
    ax = sns.stripplot(data=mrates, x="mrate", y="value", jitter=True, palette=colour, zorder=1, ax=ax, hue='Tissue', legend=False)
    ax.set_ylim(ylim)
    ax.set_xlabel("All SMs")
    ax.set_ylabel('Mutation load')
    ax.ticklabel_format(axis='y', useOffset=False, style='plain')
    yticks = ax.get_yticks()
    yticksl = yticks*10000000
    ax.set_yticks(yticks[1:-1])
    ax.set_yticklabels(yticksl[1:-1])
    ax = cleanax(ax)
    plt.tick_params(axis='x', bottom=False, labelbottom=False)

    mutsinallL1 = allsmpivot.loc[:, [b+'_L1' for b in branches]]
    mrates = pd.DataFrame(allsmpivot.loc[mutsinallL1.apply(sum, axis=1) != 7].apply(sum, axis=0)/480000000)
    mrates.columns = ['value']
    mrates['mrate'] = 'Mutation load'
    mrates['Tissue'] = [i.rsplit('_', 1)[1] for i in mrates.index.values]
    ax = fig.add_subplot(1, 2, 2)
    ax = sns.violinplot(data=mrates, x="mrate", y="value", color=colour['theme1'], linewidth=0, inner=None, ax=ax)
    plt.setp(ax.collections, alpha=.2)
    ax = sns.stripplot(data=mrates, x="mrate", y="value", jitter=True, palette=colour, zorder=1, ax=ax, hue='Tissue')
    ax.set_ylim(ylim)
    ax.set_xlabel("Filtered SMs in all L1")
    ax.set_ylabel('')
    # ax.set_ylabel('Mutation load')
    ax.ticklabel_format(axis='y', useOffset=False, style='plain')
    yticks = ax.get_yticks()
    yticksl = yticks*10000000
    ax.set_yticks(yticks[1:-1])
    ax.set_yticklabels(yticksl[1:-1])
    ax = cleanax(ax)
    ax.legend(frameon=False)
    plt.tick_params(axis='x', bottom=False, labelbottom=False)

    plt.tight_layout(pad=0)
    plt.savefig(f'{cwd}/mutation_load_differences.png', dpi=300)
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


    # <editor-fold desc="Check the presence of somatic mutations in the scRNA data">
    ## Done in python.iso_seq_analysis.get_allele_freq_at_sm_pos_plot
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


def gene_conversion_spectra():
    '''
    Read and generate summary plots from identified gene conversions
    '''
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


def axillary_stats_plots():
    """
    Functions and commands to get extra plots for manuscript
    """
    def get_fastq_readsizes(f):
        from gzip import open as gzopen
        from collections import deque
        sizes = deque()
        with gzopen(f, 'r') as fin:
            for i, line in enumerate(fin):
                if i % 4 == 1:
                    sizes.append(len(line.strip()))
        return sizes
    # END

    def raw_hifi_stats():
        """
        Get read count, read length distribution, and total number of bases
        """
        from gzip import open as gzopen
        from collections import deque
        from matplotlib import pyplot as plt
        f = '/srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4784/4784_A_run521_HIFI.fastq.gz'
        outdir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/axillary_plots/'
        sizes = get_fastq_readsizes(f)
        cnt = len(sizes)
        print(f'Number of reads: {cnt}')
        print(f'Total read lenght: {sum(sizes)}')
        plt.hist(sizes, bins=100)
        plt.ylabel("Number of reads")
        plt.xlabel("Read length")
        plt.tight_layout()
        plt.savefig(f'{outdir}/hifi_read_lenght.png', dpi=300)
        plt.close()
        return
    # END

    def get_layer_fastq_stats():
        """
            Get the sequencing read stats for the layer specific sample sequencing
        """
        from collections import defaultdict
        from glob2 import glob
        import pandas as pd
        from hometools.hometools import undict

        outdir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/axillary_plots/'
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
        # df.columns =
        df.columns = 'branch layer sequenced_reads sequenced_bases'.split()
        df['coverage'] = df['sequenced_bases']/242500000
        df.to_csv(f'{outdir}/layer_sequencing_stats.tsv', sep='\t', index=False, header=True)
        return
    # END

    def get_leaf_fastq_stats():
        """
        Get the sequencing read stats for the leaf sample sequencing
        """
        from collections import defaultdict
        from glob2 import glob
        import pandas as pd
        from hometools.hometools import undict

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
        return
    # END

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
