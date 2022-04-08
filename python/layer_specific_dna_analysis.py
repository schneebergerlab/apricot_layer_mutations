# Functions to analyse the sequencing data from the sequencing of layer-enriched dna
import matplotlib.pyplot as plt

def plot_snp_af():
    """
    This function plots the allele frequencies of the MUT_11_1 SNPs/indels in the corresponding layers
    """
    from collections import defaultdict, deque
    from matplotlib import pyplot as plt
    from matplotlib import collections as mc
    CWD = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
    LS = ['l1', 'l2', 'l3']
    MUTS = {}
    BD = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
    # Reads leaf mutations
    with open(f'{CWD}mut_11_1.mutations.regions', 'r') as fin:
        for line in fin:
            line = line.strip().split()
            MUTS[(line[0], int(line[1]))] = (line[3], line[4], int(line[5]), round(float(line[6]), 4))
    MUT_AF = defaultdict(dict)
    for k, v in MUTS.items():
        MUT_AF[k]['leaf'] = v[3]
    for l in LS:
        with open(f'{CWD}mut_11_1_{l}/read_count_at_leaf_mut_pos.txt', 'r') as fin:
            for line in fin:
                line = line.strip().split()
                if MUTS[(line[0], int(line[1]))][1] in {'A', 'C', 'G', 'T'}:
                    MUT_AF[(line[0], int(line[1]))][l] = round(int(line[BD[MUTS[(line[0], int(line[1]))][1]]])/int(line[3]), 4)
                else:
                    c = 0
                    if len(line) > 9:
                        for i in range(9, len(line), 2):
                            if line[i] == MUTS[(line[0], int(line[1]))][1]:
                                c = round(int(line[i+1])/int(line[3]), 4)
                    MUT_AF[(line[0], int(line[1]))][l] = c
    snp_lines = deque()
    indel_lines = deque()
    TYPE = ['leaf', 'l1', 'l2', 'l3']
    for k, v in MUT_AF.items():
        if MUTS[k][1] in {'A', 'C', 'G', 'T'}:
            snp_lines.append([(t, MUT_AF[k][TYPE[t]]) for t in range(4)])
        else:
            indel_lines.append([(t, MUT_AF[k][TYPE[t]]) for t in range(4)])
    # Plot SNPs and indels
    lc = mc.LineCollection(snp_lines, colors=['purple', 'lightgrey', 'purple', 'purple', 'lightgrey'] + ['purple']*4 + ['red', 'purple', 'lightgrey'], linewidths=2)
    fig = plt.figure(figsize=[5, 4])
    ax = fig.add_subplot()
    ax.set_xlim([-0.1, 3.1])
    ax.set_ylim([-0.1, 1])
    ax.add_collection(lc)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(labels=['leaf', 'layer-1', 'layer-2', 'layer-3'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Allele Frequency")
    ax.plot([], [], color='purple', label='High AF present in layers')
    ax.plot([], [], color='grey', label='Low AF absent in layers')
    ax.plot([], [], color='red', label='High AF absent in layers')
    ax.set_title("Distribution of somatic SNPs in fruits")
    plt.legend(frameon=False, bbox_to_anchor=(1, 0.95))
    plt.tight_layout()
    plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/snps_in_layers.pdf')
    plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/snps_in_layers.png', dpi=600)
    plt.close()
    lc = mc.LineCollection(indel_lines, colors='royalblue', linewidths=2)
    fig = plt.figure(figsize=[5, 4])
    ax = fig.add_subplot()
    ax.set_xlim([-0.1, 3.1])
    ax.set_ylim([-0.1, 1])
    ax.add_collection(lc)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(labels=['leaf', 'layer-1', 'layer-2', 'layer-3'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Allele Frequency")
    ax.set_title("Distribution of somatic indels in fruits")
    plt.legend(frameon=False, bbox_to_anchor=(1, 0.95))
    plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/indels_in_layers.pdf')
    plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/indels_in_layers.png', dpi=600)
    plt.close()
plot_snp_af()


def plot_snp_af_all_branch():
    """
    This function plots the allele frequencies of all SNPs/indels in the MUT_11_1 layers
    """
    from collections import defaultdict, deque
    from matplotlib import pyplot as plt
    from matplotlib import collections as mc
    CWD = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
    LS = ['l1', 'l2', 'l3']
    MUTS = {}
    BD = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
    # Reads leaf mutations
    with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mutations.regions', 'r') as fin:
        for line in fin:
            line = line.strip().split()
            MUTS[(line[0], int(line[1]))] = (line[3], line[4], int(line[5]), round(float(line[6]), 4), line[7])
    MUT_AF = defaultdict(dict)
    for k, v in MUTS.items():
        MUT_AF[k]['leaf'] = v[3]
    for l in LS:
        with open(f'{CWD}mut_11_1_{l}/read_count_at_leaf_mut_pos.all_branches.txt', 'r') as fin:
            for line in fin:
                line = line.strip().split()
                if MUTS[(line[0], int(line[1]))][1] in {'A', 'C', 'G', 'T'}:
                    MUT_AF[(line[0], int(line[1]))][l] = round(int(line[BD[MUTS[(line[0], int(line[1]))][1]]])/int(line[3]), 4)
                else:
                    c = 0
                    if len(line) > 9:
                        for i in range(9, len(line), 2):
                            if line[i] == MUTS[(line[0], int(line[1]))][1]:
                                c = round(int(line[i+1])/int(line[3]), 4)
                    MUT_AF[(line[0], int(line[1]))][l] = c
    snp_lines = deque()
    snp_sample = deque()
    indel_lines = deque()
    indel_sample = deque()
    TYPE = ['leaf', 'l1', 'l2', 'l3']
    for k, v in MUT_AF.items():
        if MUTS[k][1] in {'A', 'C', 'G', 'T'}:
            snp_lines.append([(t, MUT_AF[k][TYPE[t]]) for t in range(4)])
            snp_sample.append(MUTS[k][4])
        else:
            indel_lines.append([(t, MUT_AF[k][TYPE[t]]) for t in range(4)])
            indel_sample.append(MUTS[k][4])
    # Plot SNPs and indels
    lc = mc.LineCollection([snp_lines[i] for i in range(len(snp_sample)) if snp_sample[i] == 'MUT_11_1'], colors='purple', zorder=2, linewidths=2, label='MUT_11_1 SNPs')
    lc2 = mc.LineCollection([snp_lines[i] for i in range(len(snp_sample)) if snp_sample[i] != 'MUT_11_1'], colors='lightgreen', zorder=1, linewidths=2, label='Other SNPs')
    fig = plt.figure(figsize=[5, 4])
    ax = fig.add_subplot()
    ax.set_xlim([-0.1, 3.1])
    ax.set_ylim([-0.1, 1])
    ax.add_collection(lc)
    ax.add_collection(lc2)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(labels=['leaf', 'layer-1', 'layer-2', 'layer-3'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Allele Frequency")
    # ax.plot([], [], color='purple', label='High AF present in layers')
    # ax.plot([], [], color='grey', label='Low AF absent in layers')
    # ax.plot([], [], color='red', label='High AF absent in layers')
    ax.set_title("Distribution of somatic SNPs in fruits")
    plt.legend(frameon=False, bbox_to_anchor=(1, 0.95))
    plt.tight_layout()
    plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/snps_in_layers.all_branches.pdf')
    plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/snps_in_layers.all_branches.png', dpi=600)
    plt.close()
    lc = mc.LineCollection([indel_lines[i] for i in range(len(indel_sample)) if indel_sample[i] == 'MUT_11_1'], colors='purple', zorder=2, linewidths=2, label='MUT_11_1 indels')
    lc2 = mc.LineCollection([indel_lines[i] for i in range(len(indel_sample)) if indel_sample[i] != 'MUT_11_1'], colors='lightgreen', zorder=1, linewidths=2, label='Other indels')
    fig = plt.figure(figsize=[5, 4])
    ax = fig.add_subplot()
    ax.set_xlim([-0.1, 3.1])
    ax.set_ylim([-0.1, 1])
    ax.add_collection(lc)
    ax.add_collection(lc2)
    ax.set_xticks([0, 1, 2, 3])
    ax.set_xticklabels(labels=['leaf', 'layer-1', 'layer-2', 'layer-3'])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel("Sample")
    ax.set_ylabel("Allele Frequency")
    ax.set_title("Distribution of somatic indels in fruits")
    plt.legend(frameon=False, bbox_to_anchor=(1, 0.95))
    plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/indels_in_layers.all_branches.pdf')
    plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/indels_in_layers.all_branches.png', dpi=600)
    plt.close()
plot_snp_af_all_branch()


def layer_specific_sm_calling():
    '''
    Get somatic mutations in the three layer samples.
    Filters for: 1) positions overlapping syri SNPs/indels, 2) Positions with SM in mutiple layers
    '''

    from collections import deque, defaultdict
    from tqdm import tqdm
    ## Read snps/indels between assemblies identified by syri
    syriout = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out'
    syri_snp_list = {}
    syri_indel_list = {}
    with open(syriout, 'r') as fin:
        types = {'SNP', 'INS', 'DEL'}
        for line in fin:
            line = line.strip().split()
            if line[10] == 'SNP':
                if line[3] == 'N' or line[4] == 'N': continue
                syri_snp_list[(line[0], int(line[1]))] = (line[3], line[4])
            elif line[10] == 'INS':
                syri_indel_list[(line[0], int(line[1]))] = ('+', line[4][1:])
            elif line[10] == 'DEL':
                syri_indel_list[(line[0], int(line[1])+1)] = ('-', line[3][1:])

    ## Conditions for selecting variants:
    ## Only alt alleles with more than 5 reads are considered
    ## If a position has multiple alt_alleles => noise_pos_list AND not a somatic mutation
    ## For syri SNP/indel positions: record alt_allele_frequency AND not a somatic mutation
    ## If not above => save candidate SNP position
    indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/mut_11_1_'
    SAMPLES = ('l1', 'l2', 'l3')
    BASE_DICT = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
    BASE_DICT2 = {4: 'A', 5: 'C', 6: 'G', 7: 'T'}
    S_COV = {'l1': (50, 150) ,
             'l2': (35, 130),
            'l3': (20,180)}
    # Testing with lower read depth
    # S_COV = {'WT_1': (20, 230),
    #          'wt7': (20, 150),
    #          'wt18': (20, 150),
    #          'WT_19': (20, 240),
    #          'mut4': (20, 140),
    #          'MUT_11_1': (20, 250),
    #          'mut11_2': (20, 140),
    #          'MUT_15': (20, 230)}
    noise_pos = set()
    snp_pos = set(syri_snp_list.keys())
    indel_pos = set(syri_indel_list.keys())
    samples_sm = {}
    sample_syri_snp = {}
    sample_syri_indel = {}
    for sample in SAMPLES:
        sample_noise = deque()
        snp_alfrq = {}
        indel_alfrq = {}
        sm_pos = {}
        with open(f'{indir}{sample}/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt') as fin:
            # count = 1
            for line in tqdm(fin):
                # count += 1
                # if count == 10000:
                #     break
                line = line.strip().split()
                # Get allele frequency if position is in syri_snp
                if (line[0], int(line[1])) in snp_pos:
                    p = (line[0], int(line[1]))
                    snp_alfrq[p] = (round(int(line[BASE_DICT[syri_snp_list[p][1]]])/int(line[3]), 2), int(line[3]))
                    continue
                # Get allele frequency if position is in syri_indel
                if (line[0], int(line[1])) in indel_pos:
                    p = (line[0], int(line[1]))
                    v = syri_indel_list[p][0] + syri_indel_list[p][1]
                    for j in range(9, len(line), 2):
                        if line[j] == v:
                            indel_alfrq[p] = (round(int(line[j+1])/int(line[3]), 2), int(line[3]))
                            break
                    continue
                # Check if the position is noisy, if not then select somatic mutation
                if (line[0], int(line[1])) not in noise_pos:
                    # Check read-depth at the position
                    if not S_COV[sample][0] < int(line[3]) < S_COV[sample][1]:
                        sample_noise.append((line[0], int(line[1])))
                        continue
                    ind = [4, 5, 6, 7] + list(range(10, len(line), 2))
                    hf = 0
                    for i in ind:
                        rc = 0
                        ar = 0
                        alt = ''
                        for i in ind:
                            if int(line[i]) >= 5:
                                if BASE_DICT[line[2]] != i:
                                    if int(line[i]) <= rc:
                                        continue
                                    rc = int(line[i])
                                    ar = round(int(line[i])/int(line[3]), 2)
                                    try:
                                        alt = BASE_DICT2[i]
                                    except KeyError:
                                        alt = line[i-1]
                        if rc >= 5:
                            sm_pos[(line[0], int(line[1]))] = (line[2], alt, rc, ar)
        noise_pos = noise_pos.union(set(sample_noise))
        samples_sm[sample] = sm_pos
        sample_syri_snp[sample] = snp_alfrq
        sample_syri_indel[sample] = indel_alfrq
    # Remove noise positions from candidate lists
    for sample in SAMPLES:
        print(sample, len(samples_sm[sample]))
        for p in noise_pos:
            try:
                samples_sm[sample].pop(p)
            except KeyError:
                pass
        print(sample, len(samples_sm[sample]))

    plt.hist([v[3] for v in samples_sm['l1'].values()], bins=[i/100 for i in range(101)], histtype='step', label='L1')
    plt.hist([v[3] for v in samples_sm['l2'].values()], bins=[i/100 for i in range(101)], histtype='step', label='L2')
    plt.hist([v[3] for v in samples_sm['l3'].values()], bins=[i/100 for i in range(101)], histtype='step', label='L3')
    plt.xlabel("Allele Frequency")
    plt.ylabel("SNP/Indel Count")
    plt.legend()
    plt.tight_layout()



    # Transform data so that dict_keys are positions and dict_values are RC/AF at different levels
    snpslist = set(list(samples_sm['l1'].keys()) + list(samples_sm['l2'].keys()) + list(samples_sm['l3'].keys()))
    snpsdict = defaultdict(dict)
    for snp in snpslist:
        for sample in SAMPLES:
            try:
                snpsdict[snp][sample] = samples_sm[sample][snp]
            except KeyError:
                snpsdict[snp][sample] = (None, None, 0, 0)

    plt.close()
    afdiff = {k: max([v1[3] for v1 in v.values()]) - min([v1[3] for v1 in v.values()]) for k,v in snpsdict.items()}
    plt.hist(afdiff.values(), bins=[i/100 for i in range(101)], histtype='step', label='No RC filtering')
    afdiff = {k: max([v1[3] for v1 in v.values()]) - min([v1[3] for v1 in v.values()]) if max([v1[2] for v1 in v.values()]) > 10 else 0  for k,v in snpsdict.items()}
    plt.hist(afdiff.values(), bins=[i/100 for i in range(101)], histtype='step', label='Max RC > 10')
    plt.legend()
    plt.yscale("linear")
    plt.ylabel("# SMs")
    plt.xlabel("Max AF difference")
    plt.tight_layout()

    ## Remove positions that are present in het (0.3-0.65) in all samples
    to_pop = deque()
    for k in samples_sm['l1']:
        if not 0.3 <= samples_sm['l1'][k][3] <= 0.65:
            continue
        try:
            pop = False
            for sample in ('l2', 'l3'):
                if not 0.3 <= samples_sm[sample][k][3] <= 0.65:
                    pop = False
                    break
                else:
                    pop = True
            if pop:
                to_pop.append(k)
        except KeyError:
            pass

    ## Remove positions where afdiff (with RC filter) <0.25
    for k,v in afdiff.items():
        if v <= 0.25:
            to_pop.append(k)

    ## Many of the remaining positions are obvious mutations in all other branches. So, filter out positions that are present in other branches as well.






    for sample in SAMPLES:
        print(sample, len(samples_sm[sample]))
        for p in to_pop:
            try:
                samples_sm[sample].pop(p)
            except KeyError:
                pass
        print(sample, len(samples_sm[sample]))



    # Plot allele frequency distribution
    from matplotlib import pyplot as plt
    from matplotlib import colors as mcol
    import matplotlib

    matplotlib.rcParams['font.size'] = 8
    SAMCOLS = dict(zip(SAMPLES, [mcol.to_hex(plt.get_cmap('tab10')(i)) for i in range(len(SAMPLES))]))
    fig = plt.figure(figsize=[10, 10])
    i = 1
    for sample in SAMPLES:
        ax = fig.add_subplot(3, 2, i)
        ax.hist([v[3] for v in samples_sm[sample].values()], range=[0, 1], bins=[i/100 for i in range(101)], color=SAMCOLS[sample])
        # ax.set_yscale('log')
        ax.set_title(f'{sample} allele frequency')
        ax.axvline(x=0.25, color='black', ls=':')
        i += 1
        ax = fig.add_subplot(3, 2, i)
        ax.hist([v[2] for v in samples_sm[sample].values()], range=[0, 100], bins=[i for i in range(101)], color=SAMCOLS[sample])
        # ax.set_yscale('log')
        ax.set_title(f'{sample} read count')
        ax.axvline(x=20, color='black', ls=':')
        i += 1
    plt.tight_layout()
    plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/af_rc_after_initial_filtering.pdf')
    plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/af_rc_after_initial_filtering.png')

    # Segments positions as high_conf and low conf
    hcsm = {}        # High conf somatic mutations
    hcpos = deque()
    for sample in SAMPLES:
        sample_hcsm = {}
        for k, v in samples_sm[sample].items():
            if v[3] > 0.25 and v[2] > 20:
                sample_hcsm[k] = v
                hcpos.append(f'{k[0]}_{k[1]}')
        hcsm[sample] = sample_hcsm

    # 73 positions are selected. Manual curations still need to be done.
    # Plot hcsm using function from get_candidate_for_mutant_phenotype_funcs.py
    from get_candidate_for_mutant_phenotype_funcs import plot_selected_pos
    plot_selected_pos(hcpos, "tmp_igv.bat", "tmp/layer_specific/", M=75, HEIGHT=600, emptydir=True)
#END
layer_specific_sm_calling()

def layer_specific_gene_conversion():
#END
layer_specific_gene_conversion()


def test_allele_frequency_variation():
    '''
    Layer 3 sequencing has abnormal coverage distribution, suggesting dna degradation.
    Here, I compare the allele frequency at SNP positions in the MUT_11 leaf, layer 1, layer 2, and layer 3.
    If SNPs in L3 have similar allele frequency as other samples, then it would mean that the absurdities in
    coverage are not affecting allele frequencies and that this data can be used for mutation calling.
    '''
    from tqdm import tqdm
    from collections import deque
    from matplotlib import pyplot as plt
    import numpy as np
    import sys
    sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/software/hometools/')
    from myUsefulFunctions import density_scatter
    plt.interactive(False)

    ## Read snps between assemblies identified by syri
    syriout = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out'
    syri_snp_list = {}
    with open(syriout, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            if line[10] == 'SNP':
                if line[3] == 'N' or line[4] == 'N': continue
                syri_snp_list[(line[0], int(line[1]))] = (line[3], line[4])
    snp_pos = set(syri_snp_list.keys())
    BD = {'A': 4, 'C': 5, 'G': 6, 'T': 7}

    # Instead of re-calculating, use the pickled AFs
    # leaf_af = {}
    # with open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/MUT_11_1/bam_read_counts_b30_q10.bt2.txt", 'r') as fin:
    #     for line in tqdm(fin):
    #         line = line.strip().split()
    #         if (line[0], int(line[1])) in snp_pos:
    #             leaf_af[(line[0], int(line[1]))] = (int(line[3]) ,int(line[BD[syri_snp_list[(line[0], int(line[1]))][1]]]))
    # layer_af = deque()
    # for l in ("l1", "l2", "l3"):
    # # for l in ("l1"):
    #     ldict = {}
    #     with open(f"/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/mut_11_1_{l}/bam_read_counts_b30_q10.bt2.txt", "r") as fin:
    #         for line in tqdm(fin):
    #             line = line.strip().split()
    #             if (line[0], int(line[1])) in snp_pos:
    #                 ldict[(line[0], int(line[1]))] = (int(line[3]) ,int(line[BD[syri_snp_list[(line[0], int(line[1]))][1]]]))
    #     layer_af.append(ldict)


    import pickle
    # with open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/pickled_snp_af", "wb") as f:
    #     pickle.dump([leaf_af, layer_af], f)
    with open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/pickled_snp_af", "rb") as f:
        garb = pickle.load(f)
        leaf_af = garb[0]
        layer_af = garb[1]
    snps_in_all = set(leaf_af.keys()).intersection(set(layer_af[0].keys())).intersection(set(layer_af[1].keys())).intersection(set(layer_af[2].keys()))
    fig = plt.figure(figsize=[12,12])
    # plot Leaf vs L3
    ax1 = fig.add_subplot(3,3,1)
    ax1.set_xlim([-0.02, 1.02])
    ax1.set_ylim([-0.02, 1.02])
    x = np.array([0 if leaf_af[p][0] == 0 else leaf_af[p][1]/leaf_af[p][0] for p in snps_in_all])
    y = np.array([0 if layer_af[2][p][0] == 0 else layer_af[2][p][1]/layer_af[2][p][0] for p in snps_in_all])
    print(round(np.corrcoef(x, y)[0,1], 4))
    ax1 = density_scatter( x, y, bins = [[i/100 for i in range(0,104,4)], [i/100 for i in range(0,104,4)]], ax=ax1, fig=fig, s=0.5)
    ax1.set_xlabel("leaf")
    ax1.set_ylabel("L3")
    ax1.set_title(f"Correlation: {round(np.corrcoef(x, y)[0,1], 4)}")

    # plot Leaf vs L2
    ax2 = fig.add_subplot(3,3,4)
    ax2.set_xlim([-0.02, 1.02])
    ax2.set_ylim([-0.02, 1.02])
    x = np.array([0 if leaf_af[p][0] == 0 else leaf_af[p][1]/leaf_af[p][0] for p in snps_in_all])
    y = np.array([0 if layer_af[1][p][0] == 0 else layer_af[1][p][1]/layer_af[1][p][0] for p in snps_in_all])
    print(round(np.corrcoef(x, y)[0,1], 4))
    ax2 = density_scatter(x, y, bins = [[i/100 for i in range(0,104,4)], [i/100 for i in range(0,104,4)]], ax=ax2, fig=fig, s=0.5)
    ax2.set_xlabel("leaf")
    ax2.set_ylabel("L2")
    ax2.set_title(f"Correlation: {round(np.corrcoef(x, y)[0,1], 4)}")

    # plot Leaf vs L1
    ax3 = fig.add_subplot(3,3,7)
    ax3.set_xlim([-0.02, 1.02])
    ax3.set_ylim([-0.02, 1.02])
    x = np.array([0 if leaf_af[p][0] == 0 else leaf_af[p][1]/leaf_af[p][0] for p in snps_in_all])
    y = np.array([0 if layer_af[0][p][0] == 0 else layer_af[0][p][1]/layer_af[0][p][0] for p in snps_in_all])
    print(round(np.corrcoef(x, y)[0,1], 4))
    ax3 = density_scatter(x, y, bins = [[i/100 for i in range(0,104,4)], [i/100 for i in range(0,104,4)]], ax=ax3, fig=fig, s=0.5)
    ax3.set_xlabel("leaf")
    ax3.set_ylabel("L1")
    ax3.set_title(f"Correlation: {round(np.corrcoef(x, y)[0,1], 4)}")

    # plot L1 vs L3
    ax4 = fig.add_subplot(3, 3, 2)
    ax4.set_xlim([-0.02, 1.02])
    ax4.set_ylim([-0.02, 1.02])
    x = np.array([0 if layer_af[0][p][0] == 0 else layer_af[0][p][1]/layer_af[0][p][0] for p in snps_in_all])
    y = np.array([0 if layer_af[2][p][0] == 0 else layer_af[2][p][1]/layer_af[2][p][0] for p in snps_in_all])
    print(round(np.corrcoef(x, y)[0,1], 4))
    ax4 = density_scatter(x, y, bins = [[i/100 for i in range(0,104,4)], [i/100 for i in range(0,104,4)]], ax=ax4, fig=fig, s=0.5)
    ax4.set_xlabel("L1")
    ax4.set_ylabel("L3")
    ax4.set_title(f"Correlation: {round(np.corrcoef(x, y)[0,1], 4)}")

    # plot L1 vs L2
    ax5 = fig.add_subplot(3,3,5)
    ax5.set_xlim([-0.02, 1.02])
    ax5.set_ylim([-0.02, 1.02])
    x = np.array([0 if layer_af[0][p][0] == 0 else layer_af[0][p][1]/layer_af[0][p][0] for p in snps_in_all])
    y = np.array([0 if layer_af[1][p][0] == 0 else layer_af[1][p][1]/layer_af[1][p][0] for p in snps_in_all])
    print(round(np.corrcoef(x, y)[0,1], 4))
    ax5 = density_scatter(x, y, bins = [[i/100 for i in range(0,104,4)], [i/100 for i in range(0,104,4)]], ax=ax5, fig=fig, s=0.5)
    ax5.set_xlabel("L1")
    ax5.set_ylabel("L2")
    ax5.set_title(f"Correlation: {round(np.corrcoef(x, y)[0,1], 4)}")

    # plot L2 vs L3
    ax6 = fig.add_subplot(3, 3,3)
    ax6.set_xlim([-0.02, 1.02])
    ax6.set_ylim([-0.02, 1.02])
    x = np.array([0 if layer_af[1][p][0] == 0 else layer_af[1][p][1]/layer_af[1][p][0] for p in snps_in_all])
    y = np.array([0 if layer_af[2][p][0] == 0 else layer_af[2][p][1]/layer_af[2][p][0] for p in snps_in_all])
    print(round(np.corrcoef(x, y)[0,1], 4))
    ax6 = density_scatter(x, y, bins = [[i/100 for i in range(0,104,4)], [i/100 for i in range(0,104,4)]], ax=ax6, fig=fig, s=0.5)
    ax6.set_xlabel("2")
    ax6.set_ylabel("L3")
    ax6.set_title(f"Correlation: {round(np.corrcoef(x, y)[0,1], 4)}")
    plt.tight_layout()
    plt.savefig("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/SNP_AF_correlation.png", dpi=300)
# END
test_allele_frequency_variation()



