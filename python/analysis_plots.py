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


def mutation_spectra():
    '''
        Plots describing mutations spectra of the selected somatic mutations
    '''
    import pandas as pd
    from matplotlib import pyplot as plt
    from itertools import permutations
    import seaborn as sns
    from matplotlib.patches import Rectangle
    from collections import deque
    import numpy as np
    from pheatmap import pheatmap
    import pybedtools as bt
    import subprocess
    from scipy.cluster.hierarchy import dendrogram, linkage, optimal_leaf_ordering, leaves_list
    from scipy.sparse.csgraph import minimum_spanning_tree
    import igraph as ig
    from scipy.sparse import triu, csr_matrix
    from scipy.spatial.distance import squareform
    # from scipy.sparse import


    cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
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

    plt.rc('font', size=29.5)
    fig = plt.figure(figsize=[28, 8])
    ax = fig.add_subplot()
    sns.barplot(data=df, x='SNP type', y='Percentage of variants', hue='Layer', order=['A>G', 'T>C', 'G>A', 'C>T', 'A>C', 'T>G', 'A>T', 'T>A', 'G>T', 'C>A', 'G>C', 'C>G'], ax=ax)
    ax.legend(frameon=False, title='Layer')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout(pad=0)
    plt.savefig(f'{cwd}/snps_mutation_spectra.svg', dpi=600)
    subprocess.call(f'inkscape {cwd}/snps_mutation_spectra.svg -M {cwd}/snps_mutation_spectra.emf', shell=True)
    plt.close()

    # Heatmap comparing allele frquencies of somatic variations in leaf, L1 and L2
    ## Get readcounts from the leaf BAM files, using code in layer_specific_dna_analysis_all_samples.sh
    datafilter = data.loc[:, ['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele', 'read_count', 'allele_freq', 'Layer', 'type']].drop_duplicates()    # Filtered dataframe to be used for generating statistics
    datafilter.to_csv('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.filtered.txt', index=False, header=True, sep='\t')
    pos = set(zip(datafilter.chromosome, datafilter.position, datafilter.alt_allele))
    rc = dict()
    af = dict()
    basedict = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
    sdict = dict(zip(('WT_1', 'wt7', 'wt18', 'WT_19', 'MUT_11_1', 'mut11_2', 'mut4', 'MUT_15'), ('wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_4', 'mut_15')))
    ## Read read counts at the layer SM positions in leaves
    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
    for b in ('WT_1', 'wt7', 'wt18', 'WT_19', 'MUT_11_1', 'mut11_2', 'MUT_15', 'mut4'):
        with open(f'{cwd}/{b}.all_layer_somatic_variants.read_count.txt', 'r') as fin:
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
                    # alt = pos[line[0], int(line[1])]
                    # if alt in {'A', 'C', 'G', 'T'}:
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

    ## Read read counts at the layer SM positions in layers
    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
    for l in ('l1', 'l2', 'l3'):
        for b in ('wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15'):
            with open(f'{cwd}/{b}.{l}.all_layer_somatic_variants.read_count.txt', 'r') as fin:
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
                    # alt = pos[line[0], int(line[1])]
                    # if alt in {'A', 'C', 'G', 'T'}:
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

    rcvalues = pd.DataFrame(rc)
    afvalues = pd.DataFrame(af)
    rcvalues['var_type'] = [datafilter.loc[(datafilter.chromosome==k[0]) & (datafilter.position==k[1]) & (datafilter.alt_allele==k[2])].iloc[0].type for k in rcvalues.index.values]
    rcvalues['var_tissue'] = [datafilter.loc[(datafilter.chromosome==k[0]) & (datafilter.position==k[1]) & (datafilter.alt_allele==k[2])].iloc[0].Layer for k in rcvalues.index.values]
    rcvalues.sort_values(['var_type', 'var_tissue'], inplace=True)
    afvalues['var_type'] = [datafilter.loc[(datafilter.chromosome==k[0]) & (datafilter.position==k[1]) & (datafilter.alt_allele==k[2])].iloc[0].type for k in afvalues.index.values]
    afvalues['var_tissue'] = [datafilter.loc[(datafilter.chromosome==k[0]) & (datafilter.position==k[1]) & (datafilter.alt_allele==k[2])].iloc[0].Layer for k in afvalues.index.values]
    afvalues.sort_values(['var_type', 'var_tissue'], inplace=True)

    annocol = pd.DataFrame({
        'branch': list(rcvalues.columns[:8]) + [c.rsplit('_', maxsplit=1)[0] for c in rcvalues.columns[8:29]],
        'tissue': np.repeat(['leaf', 'L1', 'L2', 'L3'], [8, 7, 7, 7])
    })
    annorow = pd.DataFrame({
        'var_type': rcvalues.var_type,
        'var_tissue': rcvalues.var_tissue
    })
    rcvalues.to_csv(f'{cwd}/all_layer_somatic_variants.read_counts.txt', sep='\t')
    afvalues.to_csv(f'{cwd}/all_layer_somatic_variants.allele_freq.txt', sep='\t')
    annocol.to_csv(f'{cwd}/all_layer_somatic_variants.annocol.txt', sep='\t')
    annorow.to_csv(f'{cwd}/all_layer_somatic_variants.annorow.txt', sep='\t')

    # These figures did not look super great so created them using pheatmap in R (file: /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/R/presentation_plots.R)

    # fig = pheatmap(rcvalues.iloc[:, 0:29], name='ReadCount', annotation_col=annocol, annotation_row=annorow, show_rownames=False, show_colnames=False, cmap='BuGn', legend_tick_labels_styles=dict(size=10), legend_title_styles=dict(size=10), rownames_style={'fontsize':10}, legend_bar_space=1.5, legend_tick_labels={'branch': ['M1', 'M2', 'M3', 'M4', 'W1', 'W2', 'W3', 'W4']}, height=7, width=11, annotation_col_cmaps={'branch': "tab20c", 'tissue': 'tab10'}, annotation_row_cmaps={'var_type': "Accent", 'var_tissue': 'tab10'}, wspace=0.2, hspace=0.2)
    # ax = fig.get_axes()[0]
    # ax.set_xlabel('Samples')
    # ax.set_ylabel('Somatic Mutations')
    # fig.savefig(f"{cwd}/all_layer_somatic_variants.read_counts.png", dpi=600)
    # plt.close()
    # fig = pheatmap(afvalues.iloc[:, 0:29], name='AlleleFreq', annotation_col=annocol, annotation_row=annorow, show_rownames=False, show_colnames=False, cmap='BuGn', legend_tick_labels_styles=dict(size=10), legend_title_styles=dict(size=10), rownames_style={'fontsize':10}, legend_bar_space=1.5, legend_tick_labels={'branch': ['M1', 'M2', 'M3', 'M4', 'W1', 'W2', 'W3', 'W4']}, height=7, width=11, annotation_col_cmaps={'branch': "tab20c", 'tissue': 'tab10'}, annotation_row_cmaps={'var_type': "Accent", 'var_tissue': 'tab10'}, wspace=0.2, hspace=0.2)
    # ax = fig.get_axes()[0]
    # ax.set_xlabel('Samples')
    # ax.set_ylabel('Somatic Mutations')
    # fig.savefig(f"{cwd}/all_layer_somatic_variants.allele_freq.png", dpi=600)
    # plt.close()
    subprocess.call(f'inkscape {cwd}/snps_mutation_spectra.svg -M {cwd}/snps_mutation_spectra.emf', shell=True)
    subprocess.call(f'inkscape {cwd}/all_layer_somatic_variants.allele_freq.pdf -M {cwd}/all_layer_somatic_variants.allele_freq.emf', shell=True)

    # SMs overlapping genes, TEs, and intergenic, annotations
    gff = bt.BedTool('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.pasa_out.sort.protein_coding.3utr.gff3').saveas()
    te_rep = bt.BedTool('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/repeat/RepeatMasker/cur.genome.v1.fasta.ann.gff3').saveas()

    genes = gff.filter(lambda x: x[2]=='gene').saveas()
    cds = gff.filter(lambda x: x[2]=='CDS').saveas()
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
    plt.rc('font', size=14)
    fig = plt.figure(figsize=[14, 4])

    ax = fig.add_subplot(1, 2, 1)
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
    ax.set_xticklabels(['', 'Coding', 'UTR', 'Intron', 'TE/Repeat', 'Intergenic'])
    ax.legend(frameon=False, loc=2)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax = fig.add_subplot(1, 2, 2)
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
    ax.set_xticklabels(['', 'Coding', 'UTR', 'Intron', 'TE/Repeat', 'Intergenic'])
    ax.set_ylabel('Mutations per MBp')
    ax.set_xlabel('Annotation')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.tight_layout(pad=0)
    plt.savefig(f'{cwd}/all_layer_somatic_variants.annotation_overlap.svg', dpi=600)
    plt.close()
    subprocess.call(f'inkscape {cwd}/all_layer_somatic_variants.annotation_overlap.svg -M {cwd}/all_layer_somatic_variants.annotation_overlap.emf', shell=True)

    # Indel size distribution
    inddf = datafilter.loc[datafilter['type']=='Indel'].copy()
    inddf.drop_duplicates(subset=['chromosome', 'position', 'ref_allele', 'alt_allele'], inplace=True)
    inddf['size'] = [len(a)-1 if '+' in a else 1-len(a) for a in inddf.alt_allele]
    plt.hist(inddf['size'], bins=29)

    # Get maximum spanning tree and hierarchical clustering between branches (this would show
    # whether the SMs fit the morphology of the tree)
    pos = datafilter[['chromosome', 'position', 'branch', 'ref_allele', 'alt_allele']].copy()
    pos.index = pos.branch
    pos.drop(['branch'], axis=1, inplace=True)
    pos = pos.groupby(['chromosome', 'position', 'ref_allele', 'alt_allele']).filter(lambda x: len(x)<7)
    branchpos = {grp[0]: {tuple(r) for r in grp[1].itertuples(index=False)} for grp in pos.groupby(level=0)}
    branches = ['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']
    fig = plt.figure()
    ax = fig.add_subplot(1, 2, 1)

    ## Get distances for hierarchical clustering
    AM = deque()
    for i, s in enumerate(['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']):
        for j, e in enumerate(['wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15']):
            if i >= j:
                continue
            AM.append(30 - len(branchpos[s].intersection(branchpos[e])))
    AM = list(AM)
    Z = linkage(AM, method='ward')
    dendrogram(optimal_leaf_ordering(Z, AM), ax=ax)
    # ax = plt.gca()
    ax.spines[:].set_visible(False)
    ax.tick_params(left=False, labelleft=False)
    ax.set_xticks(labels=[branches[i] for i in leaves_list(optimal_leaf_ordering(Z, AM))], ticks=plt.xticks()[0], rotation=90)
    ax.set_title("Clustering of branches based on SMs")
    plt.tight_layout(pad=1)

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
    ax.set_title('Shared SMs between branches')
    ax = fig.add_subplot(2, 2, 4)
    ig.plot(T, target=ax, layout='davidson_harel', vertex_label=T.vs["name"], edge_width=np.log1p(T.es['weight']), vertex_size=0.2)
    ax.set_title('Maximum spanning tree')
    plt.tight_layout(pad=2)

    # Get AF plots of layer mutations in leaves (this would demonstrate that the
    # layer 1 mutations could not be identified in leaves, further support that
    # accurate SM identification would require homogenous cell population)
    afvalues = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.allele_freq.txt')
    afvalues.columns = ['chromosome', 'position', 'alt_allele'] + list(afvalues.columns)[3:]
    fig = plt.figure(figsize=[12, 20])
    plt.rc('font', size=10)
    for i, grp in enumerate(datafilter.groupby('branch')):
        bmut = grp[1].merge(afvalues, on=['chromosome', 'position', 'alt_allele'])
        assert(bmut.shape[0] == grp[1].shape[0])
        print(grp[0], grp[1].shape[0])
        pltdf = bmut[['chromosome', 'position', 'Layer', 'type'] + [grp[0]+s for s in ['', '_l1', '_l2', '_l3']]]
        pltdf.columns = ['chromosome', 'position', 'Layer', 'type'] + ['leaf', 'l1', 'l2', 'l3']
        pltdf = pltdf.melt(id_vars=['chromosome', 'position', 'Layer', 'type'])
        ax = fig.add_subplot(4, 2, i+1)
        sns.lineplot(pltdf, x='variable', y='value', units='position', estimator=None, hue='Layer', style='type', ax=ax)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xlabel(grp[0])
        ax.set_ylabel('Alelle freq')
        ax.set_ylim([-0.1, 1.1])
        ax.legend(loc='upper left', frameon=False)
    plt.suptitle("Allele frequency of layer specific mutations")
    plt.tight_layout(pad=4)
    plt.savefig(f'{cwd}/af_dist_layer_sm_in_leaf.pdf')
    plt.close()

    # Distribution of SMs over the genome


    # Check the presence of somatic mutations in the scRNA data





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
