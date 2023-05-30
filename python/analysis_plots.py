import argparse

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


def merge_all_SM():
    """
    Merge the list of SMs identified in leaf, layers, and all samples together
    """
    import pandas as pd
    from matplotlib import pyplot as plt
    import numpy as np
    import seaborn as sns
    from collections import deque, defaultdict, Counter
    import pybedtools as bt
    from matplotlib_venn import venn2

    # TODO: Get distibution of SMs SNPs in genomic triplets. Are CpGs enriched?

    refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
    smvars = allsmmat.reset_index()
    smvars['type'] = ['SNP' if i[0] not in '+-' else 'indel' for i in smvars['alt_allele']]
    snppos = smvars.loc[smvars.type == 'SNP']   # 1 based positions

    branches = ['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']
    l1 = set()
    l2 = set()
    for b in branches:
        pos = defaultdict(set)
        for tissue in 'leaf L1 L2'.split():
            try:
                pos[tissue] = set(snppos.loc[snppos[f'{b}_{tissue}'] == 1].index.values)
            except KeyError:
                pass
        l1 = l1.union(pos['L1'])
        l2 = l2.union(pos['L2'])
        l2 = l2.union(pos['leaf'])
    shared = l1.intersection(l2)
    l1 = l1.difference(shared)
    l2 = l2.difference(shared)


    refseq = readfasta(refcur)
    trip = deque()
    for p in snppos.itertuples(index=False):
        trip.append(refseq[p.chromosome][p.position-2: p.position+1])
    snppos['trip_reg'] = trip
    tripcnt = Counter(trip)
    tripkeys = list(tripcnt.keys())
    tripcntclean = dict()
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
    tripkeys = sorted(tripcntclean, key=lambda k: tripcntclean[k], reverse=True)
    plt.bar(tripkeys, [tripcntclean[k] for k in tripkeys])
    plt.close()


    # TODO: Correlation of SM count against branch length from the primary branching


    # TODO: Correlation of SM count against nu,ber of branching events after primary branching




    return
# END


def mutation_spectra():
    """
    Plots describing mutations spectra of the selected somatic mutations
    """
    import pandas as pd
    from matplotlib import pyplot as plt
    from itertools import permutations
    import seaborn as sns
    from collections import deque
    import numpy as np
    import pybedtools as bt
    import subprocess
    from scipy.cluster.hierarchy import dendrogram, linkage, optimal_leaf_ordering, leaves_list
    from scipy.spatial.distance import squareform
    import igraph as ig
    from scipy.sparse.csgraph import minimum_spanning_tree
    from scipy.sparse import triu, csr_matrix

    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
    bnames = ['mut_11_1', 'mut_11_2', 'wt_19', 'wt_1', 'mut_15', 'wt_18', 'wt_7']
    # Get mutation spectra (transitions vs transversions) for L1 and L2
    data = pd.read_table(f'{cwd}all_layer_somatic_variants.txt')
    ## Because reads L3 are not enriched for L3 and are mostly L2, removing them
    data.loc[data.Layer == 'L1,L3', 'Layer'] = 'L1'
    data.loc[data.Layer == 'L2,L3', 'Layer'] = 'L2'

    data['type'] = 'SNP'
    data.loc[[a[0] in '+-' for a in data.alt_allele], 'type'] = 'Indel'

    nrc = deque()
    naf = deque()
    for row in data.itertuples(index=False):
        if np.isnan(row.read_count):
            l = row.Layer
            rc = row[row._fields.index(f'{l}_RC')]
            af = row[row._fields.index(f'{l}_AF')]
            nrc.append(rc)
            naf.append(af)
        else:
            nrc.append(row.read_count)
            naf.append(row.allele_freq)
    data.read_count = nrc
    data.allele_freq = naf

    datafilter = data.loc[:, ['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele', 'read_count', 'allele_freq', 'Layer', 'type']].drop_duplicates(subset=['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele', 'Layer', 'type']).copy()    # Filtered dataframe to be used for generating statistics
    datafilter.to_csv('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.filtered.txt', index=False, header=True, sep='\t')

    # Get stats/counts
    datafilter = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.filtered.txt')

    print(f'Number of SM events: {datafilter.shape[0]}')
    print(f'Number of unique SNPs: {datafilter.loc[datafilter.type == "SNP"].drop_duplicates(subset=["chromosome", "position", "ref_allele", "alt_allele", "Layer"]).shape[0]}')
    print(f'Number of unique Indel:', datafilter.loc[datafilter.type == "Indel"] \
            .drop_duplicates(subset=["chromosome", "position", "ref_allele", "alt_allele", "Layer"]) \
            .shape[0])
    print(f'Layer specific SM counts:\n', datafilter \
          .drop_duplicates(subset=["chromosome", "position", "alt_allele"]) \
          .Layer.value_counts())
    print(f'Number of branches for each SM', Counter([grp[1].shape[0] for grp in datafilter \
                                                     .groupby(['chromosome', 'position', 'alt_allele'])]))
    print("Number of SMs per branch: ", datafilter.branch.value_counts())
    # Is the statistically significance? No
    g = datafilter.branch.value_counts()
    z = (g - np.mean(g))/np.std(g)  # Z-transform shows that No value is outside the 95% range

    print("Dimeric indels: ", {i: Counter([grp[0][2] for grp in datafilter.loc[datafilter.type == 'Indel'].groupby('chromosome position alt_allele'.split())])[i] for i in ['-AT', '-TA', '+AT', '+TA']})
    print("Dimeric indels total count: ", sum([Counter([grp[0][2] for grp in datafilter.loc[datafilter.type == 'Indel'].groupby('chromosome position alt_allele'.split())])[i] for i in ['-AT', '-TA', '+AT', '+TA']]))

    # Set colours
    colour = {('L1', 'SNP'): '#3274a1',
              ('L1', 'Indel'): '#74A6C6',
              ('L2', 'SNP'): '#e1812c',
              ('L2', 'Indel'): '#EDB07B',
              'L1': '#3274a1',
              'L2': '#e1812c'}

    # Plot for number of SMs per branch
    branches = ['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']
    branchscnt = {g[0]: g[1].branch.value_counts() for g in datafilter.groupby(['Layer', 'type'])}

    fig, ax = plt.subplots(figsize=[6, 3])
    plt.rc('font', size=15)
    y_bottom = [0]*7
    bar_width=0.5
    for l in ('L1', 'L2'):
        for t in ('SNP', 'Indel'):
            y = [branchscnt[(l, t)][b] for b in branches]
            print(y)
            # y = [cnts[x] if x in cnts else 0 for x in xticks]
            ax.bar(np.arange(1, 8), y, label=f'{l}_{t}', bottom=y_bottom, color=colour[(l, t)], width=bar_width)
            y_bottom = [y_bottom[i] + y[i] for i in range(7)]
    ax.set_xlabel("Branch")
    ax.set_ylabel("Number of SMs")
    ax.set_xticks(np.arange(1, 8))
    ax.set_xticklabels(['w1', 'w2', 'w3', 'w4', 'm1', 'm2', 'm3'])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend(bbox_to_anchor=(1.01, 1))
    plt.tight_layout(pad=0.1)
    plt.savefig(f'{cwd}/sm_counts_in_branches.png', dpi=600)
    plt.close()

    # Get transversion vs transitions plot
    unisnp = data.loc[data.type == 'SNP'].drop_duplicates(subset=['chromosome', 'position', 'ref_allele', 'alt_allele', 'Layer'])
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
    df = df.melt(id_vars=['var'])
    df.columns = ['SNP type', 'Layer', 'Percentage of variants']

    plt.rc('font', size=15)
    fig = plt.figure(figsize=[8, 3.5])
    ax = fig.add_subplot()
    sns.barplot(data=df, x='SNP type', y='Percentage of variants', hue='Layer', order=['A>G', 'T>C', 'G>A', 'C>T', 'A>C', 'T>G', 'A>T', 'T>A', 'G>T', 'C>A', 'G>C', 'C>G'], ax=ax)
    ax.legend(frameon=False, title='Layer')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.tick_params(axis='x', rotation=45)
    plt.tight_layout(pad=0)
    plt.savefig(f'{cwd}/snps_mutation_spectra.png', dpi=600)
    # subprocess.call(f'inkscape {cwd}/snps_mutation_spectra.svg -M {cwd}/snps_mutation_spectra.emf', shell=True)
    plt.close()


    # Indel size distribution
    # TODO: Add distinction between L1 and L2
    inddf = datafilter.loc[datafilter['type']=='Indel'].copy()
    inddf.drop_duplicates(subset=['chromosome', 'position', 'ref_allele', 'alt_allele'], inplace=True)
    inddf['size'] = [len(a)-1 if '+' in a else 1-len(a) for a in inddf.alt_allele]
    fig = plt.figure(figsize=[5, 3])
    ax = fig.add_subplot()
    ax.hist(inddf.loc[inddf['size'] < 0, 'size'], bins=np.linspace(-17, 0, 18), label='deletions')
    ax.hist(inddf.loc[inddf['size'] > 0, 'size'], bins=np.linspace(0, 12, 13), label='insertions')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel("Indel size")
    ax.set_ylabel("Number of indels")
    ax.legend(frameon=False)
    plt.tight_layout(pad=0)
    plt.savefig(f'{cwd}/indel_size_distribution.png', dpi=600)
    plt.close()

    # SMs overlapping genes, TEs, and intergenic, annotations
    ## TODO: Permutation test for the distribution of ShVs?
    ## TODO: Create a third plot normalised by the number of identified SM in each layer
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

    bar_width = 0.40
    xticks = ['cds', 'utr', 'intron', 'te', 'intergenic']
    plt.rc('font', size=15)
    fig = plt.figure(figsize=[8, 8])

    ax = fig.add_subplot(2, 1, 1)
    y_bottom = [0]*5
    cnts = pos.loc[(pos.Layer == 'L1') & (pos.type == 'SNP'), 'anno'].value_counts()
    y = [cnts[x] if x in cnts else 0 for x in xticks]
    ax.bar(np.arange(5) - (bar_width/2), y, label='L1_SNP', bottom=y_bottom, color='#3274a1', width=bar_width)
    y_bottom = y
    cnts = pos.loc[(pos.Layer=='L1') & (pos.type=='Indel'), 'anno'].value_counts()
    y = [cnts[x] if x in cnts else 0 for x in xticks]
    ax.bar(np.arange(5) - (bar_width/2), y, label='L1_Indel', bottom=y_bottom, color='#74A6C6', width=bar_width)
    y_bottom = [0]*5
    cnts = pos.loc[(pos.Layer=='L2') & (pos.type=='SNP'), 'anno'].value_counts()
    y = [cnts[x] if x in cnts else 0 for x in xticks]
    ax.bar(np.arange(5) + (bar_width/2), y, label='L2_SNP', bottom=y_bottom, color='#e1812c', width=bar_width)
    y_bottom = y
    cnts = pos.loc[(pos.Layer=='L2') & (pos.type=='Indel'), 'anno'].value_counts()
    y = [cnts[x] if x in cnts else 0 for x in xticks]
    ax.bar(np.arange(5) + (bar_width/2), y, label='L2_Indel', bottom=y_bottom, color='#EDB07B', width=bar_width)
    ax.set_ylabel('Mutation count')
    ax.set_xlabel('Annotation')
    ax.set_xticks(np.arange(5))
    ax.set_xticklabels(['Coding', 'UTR', 'Intron', 'TE/Repeat', 'Intergenic'])
    ax.legend(frameon=False, loc=2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax = fig.add_subplot(2, 1, 2)
    y_bottom = [0]*5
    cnts = pos.loc[(pos.Layer=='L1') & (pos.type=='SNP'), 'anno'].value_counts()
    y = [(cnts[x]*1000000)/regsize[x] if x in cnts else 0 for x in xticks]
    ax.bar(np.arange(5) - (bar_width/2), y, label='L1_SNP', bottom=y_bottom, color='#3274a1', width=bar_width)
    y_bottom = y
    cnts = pos.loc[(pos.Layer=='L1') & (pos.type=='Indel'), 'anno'].value_counts()
    y = [(cnts[x]*1000000)/regsize[x] if x in cnts else 0 for x in xticks]
    ax.bar(np.arange(5) - (bar_width/2), y, label='L1_Indel', bottom=y_bottom, color='#74A6C6', width=bar_width)
    y_bottom = [0]*5
    cnts = pos.loc[(pos.Layer=='L2') & (pos.type=='SNP'), 'anno'].value_counts()
    y = [(cnts[x]*1000000)/regsize[x] if x in cnts else 0 for x in xticks]
    ax.bar(np.arange(5) + (bar_width/2), y, label='L2_SNP', bottom=y_bottom, color='#e1812c', width=bar_width)
    y_bottom = y
    cnts = pos.loc[(pos.Layer=='L2') & (pos.type=='Indel'), 'anno'].value_counts()
    y = [(cnts[x]*1000000)/regsize[x] if x in cnts else 0 for x in xticks]
    ax.bar(np.arange(5) + (bar_width/2), y, label='L2_Indel', bottom=y_bottom, color='#EDB07B', width=bar_width)
    ax.set_xticks(np.arange(5))
    ax.set_xticklabels(['Coding', 'UTR', 'Intron', 'TE/Repeat', 'Intergenic'])
    ax.set_ylabel('Mutations per MBp')
    ax.set_xlabel('Annotation')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout(pad=0)
    plt.savefig(f'{cwd}/all_layer_somatic_variants.annotation_overlap.png', dpi=600)
    plt.close()
    subprocess.call(f'inkscape {cwd}/all_layer_somatic_variants.annotation_overlap.svg -M {cwd}/all_layer_somatic_variants.annotation_overlap.emf', shell=True)


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

    # Number of branches covered by each SM
    pos = datafilter.drop_duplicates(subset=['chromosome', 'position', 'alt_allele']).copy()
    samplecnt = {grp[0]: grp[1].shape[0] for grp in datafilter.groupby(['chromosome', 'position', 'alt_allele'])}
    pos['Branch Count'] = [samplecnt[row.chromosome, row.position, row.alt_allele] for row in pos.itertuples(index=False)]
    n_branch = pos.groupby(['Layer', 'type', 'Branch Count']).size().to_dict()
    plt.rc('font', size=15)
    fig, ax = plt.subplots(figsize=[4, 3])
    y_bottom = [0]*7
    bar_width = 0.5
    for l in ('L1', 'L2'):
        for t in ('SNP', 'Indel'):
            y = [n_branch[(l, t, i)] if (l, t, i) in n_branch else 0 for i in range(1, 8)]
            ax.bar(np.arange(1, 8), y, label=f'{l}_{t}', bottom=y_bottom, color=colour[(l, t)], width=bar_width)
            y_bottom = [y_bottom[i] + y[i] for i in range(7)]
    ax.set_xlabel("Number of Branch")
    ax.set_ylabel("Number of SMs")
    ax.set_xticks(np.arange(1, 8))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.legend(bbox_to_anchor=(1.01, 1))
    plt.tight_layout(pad=0.1)
    plt.savefig(f'{cwd}/number_sm_per_branch.png', dpi=600)
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

    # Example plot showing mutations having higher AF in L1 in all branches
    pos = ('CUR1G', 34111279, '+AT')
    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
    branches = ['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']
    readcnt = defaultdict(dict)
    for bname in branches:
        for l in ('l1', 'l2', 'l3'):
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
    fig = plt.figure(figsize=[8, 4])
    ax = fig.add_subplot()
    plt.rc('font', size=15)
    ax = pltdf.plot.bar(ax=ax)
    ax.set_ylabel("Allele Frequency")
    ax.set_title("Layer specific somatic mutation present in all branches")
    ax.legend(frameon=False)
    ax.tick_params(axis='x', rotation=0)
    ax.tick_params(axis='y', rotation=0)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout(pad=0)
    plt.savefig(f"{cwd}/layer_conserved_somatic_mutation_allele_frequency.png", dpi=600)
    plt.close()

    # Check the presence of somatic mutations in the scRNA data
    ## Done in iso_seq_analysis.py

    # Check the relationships between branching/branch lengths to SMs
    datafilter = pd.read_table(f'{cwd}/all_layer_somatic_variants.filtered.txt')
    df = datafilter.copy()
    df.branch = pd.Categorical(df.branch)
    hassm = dict()
    for grp in datafilter.groupby(['chromosome', 'position', 'alt_allele']):
        hassm[grp[0]] = dict(grp[1].branch.value_counts())
    hassm = pd.DataFrame(hassm).T.reset_index(level=[0, 1, 2])
    hassm.fillna(0, inplace=True)
    ##  Get number of SM differences between different branch points
    ### mut11.1 vs mut11.2
    sum(abs(hassm.mut_11_1 - hassm.mut_11_2)) ## 33 SMs
    ### (mut11.1 & mut11.2) vs mut15
    sum(abs((hassm.mut_11_1 * hassm.mut_11_2) - hassm.mut_15)) ## 29 SMs
    ### wt1 vs wt7
    sum(abs(hassm.wt_1 - hassm.wt_7)) ## 36 SMs
    ### (wt1 & wt7) vs wt18
    sum(abs((hassm.wt_1 * hassm.wt_7) - hassm.wt_18)) ## 7 SMs
    ### (wt1 & wt7 & wt18) vs wt19
    sum(abs((hassm.wt_1 * hassm.wt_7 * hassm.wt_18) - hassm.wt_19)) ## 48 SMs
    ### (mut11.1 & mut11.2 & mut15) vs wt19
    sum(abs((hassm.mut_11_1 * hassm.mut_11_2 * hassm.mut_15) - hassm.wt_19)) ## 44 SMs

    return1
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
        # TODO: Add this position to de novo mutation list
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
