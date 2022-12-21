# Functions to analyse the sequencing data from the sequencing of layer-enriched dna
import subprocess

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
    from collections import defaultdict, deque, Counter
    from matplotlib import pyplot as plt
    from matplotlib import collections as mc
    CWD = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
    LS = ['l1', 'l2', 'l3']
    MUTS = {}
    BD = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
    # Reads leaf mutations
    MUTS = pd.read_table("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/high_cov_mutants_sorted.all_samples.txt", skiprows=1, header=None)
    MUTS = MUTS.loc[MUTS[7] == 'Y']
    MUTS = MUTS[[0, 1, 2, 3, 4, 5, 6, 8]]
    MUTS = pd.concat([MUTS.loc[MUTS[6] == 'MUT_11_1'], MUTS.loc[MUTS[6] != 'MUT_11_1']])
    MUTS.drop_duplicates(inplace=True,subset=[0, 1])
    MUTS_dict = {}
    for row in MUTS.itertuples(index=False):
        MUTS_dict[(row[0], row[1])] = (row[2], row[3], row[4], row[5], row[6], row[7])
    MUTS = MUTS_dict
    # TODO: Finish plotting AF of all SMs in the layer samples
    MUT_AF = defaultdict(dict)
    for k, v in MUTS.items():
        MUT_AF[k]['leaf'] = v[3]
    for l in LS:
        with open(f'{CWD}mut_11_1_{l}/read_count_at_leaf_mut_pos.all_branches.updated.txt', 'r') as fin:
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
    This is the initial version set-up for MUT_11_1 only. Check layer_specific_sm_calling_all_samples for the version suitable for any sample/branch.
    '''

    from collections import deque, defaultdict
    from tqdm import tqdm
    from matplotlib import pyplot as plt
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
    # Testing with lower read-depth as this would allow selection of somatic mutations in heterozygous regions as well
    S_COV = {'l1': (20, 150),
             'l2': (20, 130),
             'l3': (20, 180)}
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
                    if not S_COV[sample][0] <= int(line[3]) <= S_COV[sample][1]:
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
        for p in noise_pos:
            try: samples_sm[sample].pop(p)
            except KeyError: pass

    plt.hist([v[3] for v in samples_sm['l1'].values()], bins=[i/100 for i in range(101)], histtype='step', label='L1')
    plt.hist([v[3] for v in samples_sm['l2'].values()], bins=[i/100 for i in range(101)], histtype='step', label='L2')
    plt.hist([v[3] for v in samples_sm['l3'].values()], bins=[i/100 for i in range(101)], histtype='step', label='L3')
    plt.xlabel("Allele Frequency")
    plt.ylabel("SNP/Indel Count")
    plt.legend()
    plt.tight_layout()

    # Transform data so that dict_keys are positions and dict_values are RC/AF at different layers
    snpslist = set(list(samples_sm['l1'].keys()) + list(samples_sm['l2'].keys()) + list(samples_sm['l3'].keys()))
    snpsdict = defaultdict(dict)
    for snp in snpslist:
        for sample in SAMPLES:
            try:
                snpsdict[snp][sample] = samples_sm[sample][snp]
            except KeyError:
                snpsdict[snp][sample] = (None, None, 0, 0)

    plt.close()
    afdiff = {k: max([v1[3] for v1 in v.values()]) - min([v1[3] for v1 in v.values()]) for k, v in snpsdict.items()}
    plt.hist(afdiff.values(), bins=[i/100 for i in range(101)], histtype='step', label='No RC filtering')
    afdiff = {k: max([v1[3] for v1 in v.values()]) - min([v1[3] for v1 in v.values()]) if max([v1[2] for v1 in v.values()]) > 10 else 0 for k, v in snpsdict.items()}
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
    for k, v in afdiff.items():
        if v <= 0.25:
            to_pop.append(k)

    for p in to_pop:
        try:
            snpsdict.pop(p)
        except KeyError:
            pass

    ## Many of the remaining positions are obvious mutations in all other branches. So, filter out positions that are present in other branches as well.
    snppos = set(snpsdict)
    snpchange = {}
    for k, v in snpsdict.items():
        for l, v1 in v.items():
            if v1[0] is None: continue
            snpchange[k] = (v1[0], v1[1])

    # write the snppos positions in file and then get unfiltered read-count at these positions in all samples
    with open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/layer_SM_candidates.txt", 'w') as fout:
        for p in snppos:
            fout.write("\t".join([p[0], str(p[1]), str(p[1])]) + "\n")
    # Get read counts by running code in layer_specific_dna_analysis.sh
    bgsnps = {}
    for sample in ('WT_1', 'wt7', 'wt18', 'WT_19', 'mut4', 'mut11_2', 'MUT_15'):
        with open(f"/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/{sample}.sm_candidate.read_count.txt", 'r') as fin:
            samplesnps = {}
            for line in fin:
                line = line.strip().split()
                if (line[0], int(line[1])) in snppos:
                    rc, af = 0, 0
                    if snpchange[(line[0], int(line[1]))][1] in {'A', 'C', 'G', 'T'}:
                        rc = int(line[BASE_DICT[snpchange[(line[0], int(line[1]))][1]]])
                        af = round(int(line[BASE_DICT[snpchange[(line[0], int(line[1]))][1]]])/int(line[3]), 2)
                    else:
                        if len(line) > 9:
                            for i in range(9, len(line), 2):
                                if line[i] == snpchange[(line[0], int(line[1]))][1]:
                                    rc = int(line[i+1])
                                    af = round(int(line[i+1])/int(line[3]), 2)
                    samplesnps[(line[0], int(line[1]))] = (rc, af)
            bgsnps[sample] = samplesnps
    # Filter positions supported by at least 5 reads in all samples
    bgsnpscnt = Counter([k for v in bgsnps.values() for k, v1 in v.items() if v1[0] >= 5])
    bgsnpsuni = set([k for k, v in bgsnpscnt.items() if v == 7])
    snppos_bgfilt = snppos - bgsnpsuni
    # len(snppos_bgfilt): 297

    # Filter positions supported by at least 5 reads in the illumina libraries (this removes library specific noisy reads)
    bgsnpscnt = Counter([k for s in ('wt7', 'wt18', 'mut4', 'mut11_2') for k, v1 in bgsnps[s].items() if v1[0] >= 5])
    bgsnpsuni = set([k for k, v in bgsnpscnt.items() if v == 4])
    snppos_bgfilt = snppos_bgfilt - bgsnpsuni
    # len(snppos_bgfilt): 197

    # Select high conf candidates
    hcsm = {}        # High conf somatic mutations
    for p in snppos_bgfilt:
        for k, v in snpsdict[p].items():
            if v[3] > 0.25 and v[2] > 20:
                hcsm[p] = snpsdict[p]
                break
    # 52 hcsm of which 51 are usable
    with open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/high_conf_layer_specific_somatic_mutations.txt", 'w') as fout:
        for k in sorted(hcsm):
            v = hcsm[k]
            if len(set([v1[1] for v1 in v.values() if v1[1] is not None])) > 1: continue
            ref, alt = None, None
            for l in ['l1', 'l2', 'l3']:
                if v[l][0] is not None:
                    ref, alt = v[l][:2]
                    break
            vals = deque()
            for l in ['l1', 'l2', 'l3']:
                vals.extend(v[l][2:])
            fout.write("\t".join(list(map(str, list(k) +  [ref, alt] + list(vals)))) + "\n")
    # Do manual curation of the candidate somatic mutation list
    # 51 positions tested, 41 selected as somatic mutations (11 of these were present in leaf, 30 are new).
    # Out of 41 SM, 25 are in L1, 13 are shared by L2-L3, 1 in L1-L3, and 2 in L3.
    # 34 SMs are in diploid regions, 4 in haploid regions (NOTAL/HDR regions in syri output), and 3 appear to be in duplicated regions
#END
layer_specific_sm_calling()


def layer_specific_sm_calling_all_samples(cwd, bname, scov, nc=1):
    '''
    This is a more generalised version of the layer_specific_sm_calling() and it created so that it can be used on any branch.

    Get somatic mutations in the three layer samples.
    Filters for: 1) positions overlapping syri SNPs/indels, 2) Positions with SM in mutiple layers

    Inputs:
    1) cwd = path to folder containing different layers from a branch (e.g. /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/wt_1)
    2) branch name = e.g. 'wt_1'
    3) scov = Dictionary containing optimal read-depth ranges (e.g. {'l1': (20, 150), 'l2': (20, 130), 'l3': (20, 180)} )

    Run on DELL-NODES as it will not work on local computer due to memory limitations.

    Coverage dicts for samples:
    MUT_11_1 : {'l1': (20, 160), 'l2': (20, 150), 'l3': (20, 180)}
    '''

    from collections import deque, defaultdict, Counter
    from matplotlib import pyplot as plt
    from pandas import read_table
    from subprocess import Popen, PIPE
    import numpy as np

    ## Read snps/indels between assemblies identified by syri
    syri_snp_list, syri_indel_list = getsyrivarlist()

    ## Conditions for selecting variants:
    ## Only alt alleles with more than 5 reads are considered
    ## If a position has multiple alt_alleles => noise_pos_list AND not a somatic mutation
    ## For syri SNP/indel positions: record alt_allele_frequency AND not a somatic mutation
    ## If not above => save candidate SNP position
    indir = f'{cwd}/{bname}_' #/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/mut_11_1_'
    noise_pos = set()
    # snp_pos = set(syri_snp_list.keys())
    # indel_pos = set(syri_indel_list.keys())
    samples_sm = {}
    sample_syri_snp = {}
    sample_syri_indel = {}
    SAMPLES = ('l1', 'l2', 'l3')
    for sample in SAMPLES:
        filename = f'{indir}{sample}/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt'
        sample_noise, snp_alfrq, indel_alfrq, sm_pos = readfilteredbamreadcount(filename, scov[sample][0], scov[sample][1], noise_pos, syri_snp_list, syri_indel_list)

        noise_pos = noise_pos.union(set(sample_noise))
        samples_sm[sample] = sm_pos
        sample_syri_snp[sample] = snp_alfrq
        sample_syri_indel[sample] = indel_alfrq

    # Remove noise positions from candidate lists
    for sample in SAMPLES:
        for p in noise_pos:
            try:
                samples_sm[sample].pop(p)
            except KeyError:
                pass

    # TODO: Check and save this figure
    plt.hist([v[3] for v in samples_sm['l1'].values()], bins=[i/100 for i in range(101)], histtype='step', label='L1')
    plt.hist([v[3] for v in samples_sm['l2'].values()], bins=[i/100 for i in range(101)], histtype='step', label='L2')
    plt.hist([v[3] for v in samples_sm['l3'].values()], bins=[i/100 for i in range(101)], histtype='step', label='L3')
    plt.xlabel("Allele Frequency")
    plt.ylabel("SNP/Indel Count")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'{cwd}/sm_allele_frequency_layers.pdf')
    plt.close()

    # Transform data so that dict_keys are positions and dict_values are RC/AF at different layers
    snpslist = set(list(samples_sm['l1'].keys()) + list(samples_sm['l2'].keys()) + list(samples_sm['l3'].keys()))
    snpsdict = defaultdict(dict)
    for snp in snpslist:
        for sample in SAMPLES:
            try:
                snpsdict[snp][sample] = samples_sm[sample][snp]
            except KeyError:
                snpsdict[snp][sample] = (None, None, 0, 0)

    # TODO: Check and save this figure
    afdiff = {k: max([v1[3] for v1 in v.values()]) - min([v1[3] for v1 in v.values()]) for k, v in snpsdict.items()}
    plt.hist(afdiff.values(), bins=[i/100 for i in range(101)], alpha=0.5, label='No RC filtering')
    afdiff = {k: max([v1[3] for v1 in v.values()]) - min([v1[3] for v1 in v.values()]) if max([v1[2] for v1 in v.values()]) > 10 else 0 for k, v in snpsdict.items()}
    afcut = np.quantile(list(afdiff.values()), 0.99)

    plt.hist(afdiff.values(), bins=[i/100 for i in range(101)], alpha=0.5, label='Max RC > 10')
    plt.axvline(afcut)
    plt.legend()
    plt.yscale("linear")
    plt.ylabel("# SMs")
    plt.xlabel("Max AF difference")
    plt.tight_layout()
    plt.savefig(f'{cwd}/max_allele_frequency_between_layers.pdf')
    plt.close()

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
    for k, v in afdiff.items():
        # if v <= 0.25:
        if v <= afcut: # Testing with a smalled afdiff cutoff
            to_pop.append(k)

    for p in to_pop:
        try:
            snpsdict.pop(p)
        except KeyError:
            pass

    ## Many of the remaining positions are obvious mutations in all other branches. So, filter out positions that are present in other branches as well.
    snppos = set(snpsdict)
    snpchange = {}
    for k, v in snpsdict.items():
        for l, v1 in v.items():
            if v1[0] is None: continue
            snpchange[k] = (v1[0], v1[1])

    # write the snppos positions in file and then get unfiltered read-count at these positions in all samples
    with open(f"{cwd}/layer_SM_candidates.txt", 'w') as fout:
        for p in snppos:
            fout.write("\t".join([p[0], str(p[1]), str(p[1])]) + "\n")

    # Sort the candidates
    df = read_table(f"{cwd}/layer_SM_candidates.txt", header=None)
    df.sort_values([0, 1, 2], inplace=True)
    df.to_csv(f"{cwd}/layer_SM_candidates.txt", index=False, header=False, sep='\t')

    # Get read counts at candidate positions in the leaf samples
    leafdir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
    refgen = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
    procs = deque()
    for leaf in ('WT_1', 'WT_19', 'MUT_15', 'MUT_11_1'):
        command = f"/srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/python /srv/biodata/dep_mercier/grp_schneeberger/software/hometools/myUsefulFunctions.py pbamrc -n {nc} -b 0 -q 0 -w 0 -I -f {refgen} -l {cwd}layer_SM_candidates.txt {leafdir}/{leaf}/{leaf}.sorted.bt2.bam {cwd}{leaf}.sm_candidate.read_count.txt"
        p = Popen(command.split(), stderr=PIPE, stdout=PIPE, text=True)
        procs.append(p)
    for leaf in ('wt7', 'wt18', 'mut4', 'mut11_2'):
        command = f"/srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/python /srv/biodata/dep_mercier/grp_schneeberger/software/hometools/myUsefulFunctions.py pbamrc -n {nc} -b 0 -q 0 -w 0 -I -f {refgen} -l {cwd}layer_SM_candidates.txt {leafdir}/{leaf}/{leaf}.deduped.bam {cwd}{leaf}.sm_candidate.read_count.txt"
        p = Popen(command.split(), stderr=PIPE, stdout=PIPE, text=True)
        procs.append(p)
    for p in procs:
        o = p.communicate()
        if o[1] != '': print(o)

    # Read background read counts at SM candidate positions
    sample_ids = dict(zip(('wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_4', 'mut_15'), ('WT_1', 'wt7', 'wt18', 'WT_19', 'MUT_11_1', 'mut11_2', 'mut4', 'MUT_15')))
    sample_ids.pop(bname)
    bgsnps = {}
    # for sample in ('WT_1', 'wt7', 'wt18', 'WT_19', 'mut4', 'mut11_2', 'MUT_15'):
    for sample in sample_ids.values():
        with open(f"{cwd}/{sample}.sm_candidate.read_count.txt", 'r') as fin:
            samplesnps = {}
            for line in fin:
                line = line.strip().split()
                if (line[0], int(line[1])) in snppos:
                    rc, af = 0, 0
                    if snpchange[(line[0], int(line[1]))][1] in {'A', 'C', 'G', 'T'}:
                        rc = int(line[BASE_DICT[snpchange[(line[0], int(line[1]))][1]]])
                        af = round(int(line[BASE_DICT[snpchange[(line[0], int(line[1]))][1]]])/int(line[3]), 2)
                    else:
                        if len(line) > 9:
                            for i in range(9, len(line), 2):
                                if line[i] == snpchange[(line[0], int(line[1]))][1]:
                                    rc = int(line[i+1])
                                    af = round(int(line[i+1])/int(line[3]), 2)
                    samplesnps[(line[0], int(line[1]))] = (rc, af)
            bgsnps[sample] = samplesnps
    # Filter positions supported by at least 5 reads in all samples
    bgsnpscnt = Counter([k for v in bgsnps.values() for k, v1 in v.items() if v1[0] >= 5])
    bgsnpsuni = set([k for k, v in bgsnpscnt.items() if v == 7])
    snppos_bgfilt = snppos - bgsnpsuni
    # len(snppos_bgfilt)

    # Filter positions supported by at least 5 reads in the illumina libraries (this removes library specific noisy reads)
    slist = set([v for v in sample_ids.values() if v in {'WT_1', 'WT_19', 'MUT_11_1', 'MUT_15'}])
    bgsnpscnt = Counter([k for s in slist for k, v1 in bgsnps[s].items() if v1[0] >= 5])
    bgsnpsuni = set([k for k, v in bgsnpscnt.items() if v == len(slist)])
    snppos_bgfilt = snppos_bgfilt - bgsnpsuni
    # len(snppos_bgfilt)

    # Filter positions supported by at least 5 reads in the illumina libraries (this removes library specific noisy reads)
    slist = set([v for v in sample_ids.values() if v in {'wt7', 'wt18', 'mut4', 'mut11_2'}])
    bgsnpscnt = Counter([k for s in slist for k, v1 in bgsnps[s].items() if v1[0] >= 5])
    bgsnpsuni = set([k for k, v in bgsnpscnt.items() if v == len(slist)])
    snppos_bgfilt = snppos_bgfilt - bgsnpsuni
    # len(snppos_bgfilt)

    # Select high conf candidates
    hcsm = {}        # High conf somatic mutations
    for p in snppos_bgfilt:
        for k, v in snpsdict[p].items():
            # if v[3] > 0.25 and v[2] > 20:
            if v[3] > afcut and v[2] > 20:
                hcsm[p] = snpsdict[p]
                break
    with open(f"{cwd}/high_conf_layer_specific_somatic_mutations.txt", 'w') as fout:
        for k in sorted(hcsm):
            v = hcsm[k]
            if len(set([v1[1] for v1 in v.values() if v1[1] is not None])) > 1: continue
            ref, alt = None, None
            for l in ['l1', 'l2', 'l3']:
                if v[l][0] is not None:
                    ref, alt = v[l][:2]
                    break
            vals = deque()
            for l in ['l1', 'l2', 'l3']:
                vals.extend(v[l][2:])
            fout.write("\t".join(list(map(str, list(k) + [ref, alt] + list(vals)))) + "\n")
    # Do manual curation of the candidate somatic mutation list
#END

# wt_1
layer_specific_sm_calling_all_samples('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/wt_1/', 'wt_1', {'l1': (40, 240), 'l2': (20, 180), 'l3': (30, 240)}, nc=6)

# wt_7
layer_specific_sm_calling_all_samples('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/wt_7/', 'wt_7', {'l1': (30, 180), 'l2': (30, 180), 'l3': (30, 180)}, nc=6)

# wt_18
layer_specific_sm_calling_all_samples('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/wt_18/', 'wt_18', {'l1': (40, 200), 'l2': (30, 180), 'l3': (30, 180)}, nc=6)

# wt_19
layer_specific_sm_calling_all_samples('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/wt_19/', 'wt_19', {'l1': (40, 220), 'l2': (40, 220), 'l3': (40, 220)}, nc=6)

# mut_11_1
layer_specific_sm_calling_all_samples('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/mut_11_1/', 'mut_11_1', {'l1': (20, 160), 'l2': (20, 150), 'l3': (20, 180)}, nc=6)

# mut_11_2
layer_specific_sm_calling_all_samples('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/mut_11_2/', 'mut_11_2', {'l1': (30, 200), 'l2': (20, 180), 'l3': (30, 200)}, nc=6)

# mut_15
layer_specific_sm_calling_all_samples('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/mut_15/', 'mut_15', {'l1': (40, 220), 'l2': (30, 220), 'l3': (40, 220)}, nc=6)


def layer_conserved_variants():
    '''
    Some layer-specific variants are conserved in multi-branches. These could not
    be identified in the layer_specific_sm_calling_all_samples because for such
    positions the leaf samples have ALT reads too, making them look like noise.

    Here, I try to call variations that are present in a specific layer in multiple
    branches. For this, I would skip filtering based on leaf and try to compare
    AFs in the layers directly using t-test.
    '''

    from collections import deque, defaultdict, Counter
    from matplotlib import pyplot as plt
    from pandas import read_table
    from subprocess import Popen, PIPE
    import numpy as np
    from scipy.stats import ttest_rel

    ## Read snps/indels between assemblies identified by syri
    syri_snp_list, syri_indel_list = getsyrivarlist()

    ## Conditions for selecting variants:
    ## Only alt alleles with more than 5 reads are considered
    ## If a position has multiple alt_alleles => noise_pos_list AND not a somatic mutation
    ## For syri SNP/indel positions: record alt_allele_frequency AND not a somatic mutation
    ## If not above => save candidate SNP position
    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/'
    branches = ('wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15')
    scovs = {'wt_1':  {'l1': (40, 240), 'l2': (20, 180), 'l3': (30, 240)},
             'wt_7':  {'l1': (30, 180), 'l2': (30, 180), 'l3': (30, 180)},
             'wt_18': {'l1': (40, 200), 'l2': (30, 180), 'l3': (30, 180)},
             'wt_19': {'l1': (40, 220), 'l2': (40, 220), 'l3': (40, 220)},
             'mut_11_1': {'l1': (20, 160), 'l2': (20, 150), 'l3': (20, 180)},
             'mut_11_2': {'l1': (30, 200), 'l2': (20, 180), 'l3': (30, 200)},
             'mut_15': {'l1': (40, 220), 'l2': (30, 220), 'l3': (40, 220)}}


    branch_vars = {}
    for bname in branches:
        indir = f'{cwd}/{bname}/{bname}_'
        noise_pos = set()
        samples_sm = {}
        sample_syri_snp = {}
        sample_syri_indel = {}
        SAMPLES = ('l1', 'l2', 'l3')
        for sample in SAMPLES:
            filename = f'{indir}{sample}/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt'
            sample_noise, snp_alfrq, indel_alfrq, sm_pos = readfilteredbamreadcount(filename, scovs[bname][sample][0], scovs[bname][sample][1], noise_pos, syri_snp_list, syri_indel_list)
            noise_pos = noise_pos.union(set(sample_noise))
            samples_sm[sample] = sm_pos
            sample_syri_snp[sample] = snp_alfrq
            sample_syri_indel[sample] = indel_alfrq

        # Remove noise positions from candidate lists
        for sample in SAMPLES:
            for p in noise_pos:
                try:
                    samples_sm[sample].pop(p)
                except KeyError:
                    pass

        # Transform data so that dict_keys are positions and dict_values are RC/AF at different layers
        snpslist = set(list(samples_sm['l1'].keys()) + list(samples_sm['l2'].keys()) + list(samples_sm['l3'].keys()))
        snpsdict = defaultdict(dict)
        for snp in snpslist:
            for sample in SAMPLES:
                try:
                    snpsdict[snp][sample] = samples_sm[sample][snp]
                except KeyError:
                    snpsdict[snp][sample] = (None, None, 0, 0)

        # # TODO: Check and save this figure
        # # afdiff = {k: max([v1[3] for v1 in v.values()]) - min([v1[3] for v1 in v.values()]) for k, v in snpsdict.items()}
        # afdiff = {k: max([v1[3] for v1 in v.values()]) - min([v1[3] for v1 in v.values()]) if max([v1[2] for v1 in v.values()]) > 10 else 0 for k, v in snpsdict.items()}
        # afcut = np.quantile(list(afdiff.values()), 0.99)
        #
        # ## Remove positions that are present in het (0.3-0.65) in all samples
        # to_pop = deque()
        # for k in samples_sm['l1']:
        #     if not 0.3 <= samples_sm['l1'][k][3] <= 0.65:
        #         continue
        #     try:
        #         pop = False
        #         for sample in ('l2', 'l3'):
        #             if not 0.3 <= samples_sm[sample][k][3] <= 0.65:
        #                 pop = False
        #                 break
        #             else:
        #                 pop = True
        #         if pop:
        #             to_pop.append(k)
        #     except KeyError:
        #         pass
        #
        # ## Remove positions where afdiff (with RC filter) <0.25
        # for k, v in afdiff.items():
        #     if v <= afcut: # Testing with a smalled afdiff cutoff
        #         to_pop.append(k)
        # for p in to_pop:
        #     try:
        #         snpsdict.pop(p)
        #     except KeyError:
        #         pass
        branch_vars[bname] = snpsdict

    pos = set()
    for k, v in branch_vars.items():
        pos = pos.union(set(v.keys()))

    posdict = dict()
    for p in pos:
        posdict[p] = {}
        for bname in branches:
            v = branch_vars[bname][p]
            if len(v) != 0:
                posdict[p][bname] = v

    ## Filetering posdict
    to_pop = deque()
    ### Remove positions which are present in only 1 branch
    for k, v in posdict.items():
        if len(v) == 1:
            to_pop.append(k)

    ### Remove positions which are not supported by atleast 20 reads in any l1/l2
    for k, v in tqdm(posdict.items()):
        if any([v1['l1'][2]>20 or v1['l2'][2]>20 for v1 in v.values()]): continue
        to_pop.append(k)

    for p in to_pop:
        try:
            posdict.pop(p)
        except KeyError: pass

    ### Get log2(fold-change) values
    fc = dict()
    for k, v in posdict.items():
        a = np.mean([v1['l1'][3] for v1 in v.values()])
        b = np.mean([v1['l2'][3] for v1 in v.values()])
        if b != 0:
            if a != 0:
                f = np.log2(a/b)
            else:
                f = -10
        else:
            f = 10
        fc[k] = f

    plt.hist(fc.values(), bins=100)
    plt.yscale('log')
    plt.xlabel('Log2(FC)')
    plt.ylabel('Positions Count')
    plt.axvline(2, color='black')
    plt.axvline(-2, color='black')
    ## Select positions with log2(FC) > 2 as candidate multi-branch layer-specific variations
    selected = [k for k, v in fc.items() if abs(v) > 2]

    with open(f'{cwd}/conserved_layer_specific_candidate.tsv', 'w') as fout:
        for s in sorted(selected):
            fout.write(f'{s[0]}\t{s[1]}\t')
            st = True
            for k, v in posdict[s].items():
                if st:
                    fout.write(f'{k}\t')
                    st = False
                else:
                    fout.write(f'\t\t{k}\t')
                st2 = True
                for k1, v1 in v.items():
                    if st2:
                        fout.write(f'{k1}\t{v1[0]}\t{v1[1]}\t{v1[2]}\t{v1[3]}\n')
                        st2=False
                    else:
                        fout.write(f'\t\t\t{k1}\t{v1[0]}\t{v1[1]}\t{v1[2]}\t{v1[3]}\n')
            fout.write('\n')
    ## 64 positions selected as candidate mutations. Next I test them manually.







    # posstd = dict()
    # for k, v in posdict.items():
    #     values = [v1['l1'][3] for v1 in v.values()] +  [v1['l2'][3] for v1 in v.values()]
    #     posstd[k] = np.var(values)
    #     break
    # cutoff = np.quantile(list(posstd.values()), 0.99)
    # to_pop = deque()
    # for k, v in posstd.items():
    #     if v < cutoff:
    #         to_pop.append(k)
    #
    # for p in to_pop:
    #     try:
    #         posdict.pop(p)
    #     except KeyError: pass
    #
    # figure = plt.figure(figsize=[8, 6])
    # for i in range(2, 8):
    #     pvalues = dict()
    #     pos = [k for k, v in posdict.items() if len(v) == i]
    #     for p in pos:
    #         v = posdict[p]
    #         l1s = [v1['l1'][3] for v1 in v.values()]
    #         l2s = [v1['l2'][3] for v1 in v.values()]
    #         pvalues[p] = ttest_rel(l1s, l2s)
    #     a = np.array([v.pvalue for v in pvalues.values()])
    #     np.nan_to_num(a, copy=False, nan=1)
    #     b = p_adjust(a, 'bh')
    #     selected = np.where(np.array(b) < 0.05)[0]
    #     print(i, len(a), len(selected))
    #     hc = deque()
    #     for s in selected:
    #         l = 'l1' if sum([v['l1'][3] for v in posdict[pos[s]].values()]) > sum([v['l2'][3] for v in posdict[pos[s]].values()]) else 'l2'
    #         if any([v[l][2] > 20 for v in posdict[pos[s]].values()]):
    #             # print(pos[s], posdict[pos[s]])
    #             # break
    #             hc.append(s)
    #
    #         # print(pos[s], posdict[pos[s]], sum([v['l1'][2] for v in posdict[pos[s]].values()])/i, sum([v['l2'][2] for v in posdict[pos[s]].values()])/i)
    #         # print()
    #     print(f'i: {i}; Count: {len(a)}; Count Selected: {len(selected)}; Count high read count: {len(hc)}')
    # break




        # garb = [('CUR1G', 10626876),
        #         ('CUR1G', 2854033),
        #         ('CUR4G', 6760960),
        #         ('CUR1G', 34111279),
        #         ('CUR4G', 10315321),
        #         ('CUR5G', 10982419),
        #         ('CUR7G', 4335333),
        #         ('CUR3G', 16974418),
        #         ('CUR1G', 2981550),
        #         ('CUR5G', 16975416)]




# END



def layer_3_variant_calling(cwd, bname, scov, nc=1):
    '''
    Layer 3 could not be so well enriched, i.e. there are quite some L2 cells in L3
    sample. This mean that calling variations present in only L3 is more challenging.
    Here, try to get all variations where the variant in L3 but absent in L2 and
    then try to call L3 specific variants from there.
    '''
    from collections import defaultdict, deque, Counter
    import numpy as np
    import pybedtools as bd
    from subprocess import Popen, PIPE

    BASE_DICT = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
    BASE_DICT2 = {4: 'A', 5: 'C', 6: 'G', 7: 'T'}

    ## Read snps/indels between assemblies identified by syri
    syri_snp_list, syri_indel_list = getsyrivarlist()

    ## Conditions for selecting variants:
    ## Only alt alleles with more than 5 reads are considered
    ## If a position has multiple alt_alleles => noise_pos_list AND not a somatic mutation
    ## For syri SNP/indel positions: record alt_allele_frequency AND not a somatic mutation
    ## If not above => save candidate SNP position
    indir = f'{cwd}/{bname}_' #/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/mut_11_1/mut_11_1_'
    noise_pos = set()
    samples_sm = {}
    sample_syri_snp = {}
    sample_syri_indel = {}
    SAMPLES = ('l1', 'l2', 'l3')
    for sample in SAMPLES:
        filename = f'{indir}{sample}/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt'
        sample_noise, snp_alfrq, indel_alfrq, sm_pos = readfilteredbamreadcount(filename, scov[sample][0], scov[sample][1], noise_pos, syri_snp_list, syri_indel_list)
        noise_pos = noise_pos.union(set(sample_noise))
        samples_sm[sample] = sm_pos
        sample_syri_snp[sample] = snp_alfrq
        sample_syri_indel[sample] = indel_alfrq

    # Filter out noise positions
    for sample in SAMPLES:
        for p in noise_pos:
            try: samples_sm[sample].pop(p)
            except KeyError: pass

    # Transform data so that dict_keys are positions and dict_values are RC/AF at different layers
    snpslist = set(list(samples_sm['l1'].keys()) + list(samples_sm['l2'].keys()) + list(samples_sm['l3'].keys()))
    snpsdict = defaultdict(dict)
    for snp in snpslist:
        for sample in SAMPLES:
            try: snpsdict[snp][sample] = samples_sm[sample][snp]
            except KeyError: snpsdict[snp][sample] = (None, None, 0, 0)

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

    for p in to_pop:
        try:
            snpsdict.pop(p)
        except KeyError:
            pass

    afdiff = {k: v['l3'][3]-v['l2'][3] for k, v in snpsdict.items()}
    afcut = np.quantile(list(afdiff.values()), 0.99)

    highl3af = {k: v for k, v in snpsdict.items() if afdiff[k] > afcut and v['l3'][2]>=10 and v['l2'][2]<=1}

    goodcov = bd.BedTool('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/mut_11_1/mut_11_1_l3/l3.depth_30_300.bed')
    highl3bed = bd.BedTool('\n'.join(['\t'.join([i[0], str(i[1]-1), str(i[1])]) for i in highl3af]), from_string=True)
    garb = highl3bed.intersect(goodcov)

    with open(f"{cwd}/layer3_SM_candidates.txt", 'w') as fout:
        for g in garb:
            fout.write(f'{g[0]}\t{g[2]}\t{g[2]}\n')

    # Get read counts at candidate L3 positions in other branches with mapping quality 0
    refgen = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
    procs = deque()
    for l in ('l1', 'l2'):
        command = f"/srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/hometools pbamrc -n {nc} -b 0 -q 0 -w 0 -I -f {refgen} -l {cwd}/layer3_SM_candidates.txt {cwd}/{bname}_{l}/{l}.deduped.bam {cwd}/{l}.for_l3.sm_candidate.read_count.txt"
        p = Popen(command.split(), stderr=PIPE, stdout=PIPE, text=True)
        procs.append(p)
    for p in procs:
        o = p.communicate()
        if o[1] != '': print(o)

    # Read background read counts at SM candidate positions
    bgsnps = {}
    snpchange = {}
    for k, v in highl3af.items():
        for l, v1 in v.items():
            if v1[0] is None: continue
            snpchange[k] = (v1[0], v1[1])

    highl3pos = set([(g[0], int(g[2])) for g in garb])
    for l in ('l1', 'l2'):
        with open(f"{cwd}/{l}.for_l3.sm_candidate.read_count.txt", 'r') as fin:
            samplesnps = {}
            for line in fin:
                line = line.strip().split()
                if (line[0], int(line[1])) in highl3pos:
                    rc, af = 0, 0
                    if snpchange[(line[0], int(line[1]))][1] in {'A', 'C', 'G', 'T'}:
                        rc = int(line[BASE_DICT[snpchange[(line[0], int(line[1]))][1]]])
                        af = round(int(line[BASE_DICT[snpchange[(line[0], int(line[1]))][1]]])/int(line[3]), 2)
                    else:
                        if len(line) > 9:
                            for i in range(9, len(line), 2):
                                if line[i] == snpchange[(line[0], int(line[1]))][1]:
                                    rc = int(line[i+1])
                                    af = round(int(line[i+1])/int(line[3]), 2)
                    samplesnps[(line[0], int(line[1]))] = (rc, af)
            bgsnps[l] = samplesnps

    # Filter positions supported by at least 5 reads in both l1 and l2 samples
    bgsnpscnt = Counter([k for v in bgsnps.values() for k, v1 in v.items() if v1[0] > 1])
    bgsnpsuni = set([k for k, v in bgsnpscnt.items() if v==2])
    highl3pos_bgfilt = highl3pos - bgsnpsuni
    # len(highl3pos_bgfilt)

    # print(bname, len(highl3pos_bgfilt), highl3pos_bgfilt)
    with open(f"{cwd}/l3_specific_somatic_mutations.txt", 'w') as fout:
        for k in sorted(highl3pos_bgfilt):
            v = highl3af[k]
            if len(set([v1[1] for v1 in v.values() if v1[1] is not None])) > 1: continue
            ref, alt = None, None
            for l in ['l1', 'l2', 'l3']:
                if v[l][0] is not None:
                    ref, alt = v[l][:2]
                    break
            vals = deque()
            for l in ['l1', 'l2', 'l3']:
                vals.extend(v[l][2:])
            fout.write("\t".join(list(map(str, list(k) + [ref, alt] + list(vals)))) + "\n")

    '''
    # I checked these candidates in all samples manually. Other than 1 position in MUT_11_1, none of the candidates could be considered as L3 specific mutation. Possible reasons: 
    1) There is not much enrichment in L3 DNA, i.e. there are many L2 reads in L3 sample as a result we cannot get L3 specific variants that have enough reads/allele-frequency to be easily called. We know that L3 layer is very thin in fruit, so this is quite possible
    2) There are not many mutations happening in L3 (but this is very risky to say)
    '''
# END
# wt_1
layer_3_variant_calling('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/wt_1/', 'wt_1', {'l1': (40, 240), 'l2': (20, 180), 'l3': (30, 240)}, nc=6)

# wt_7
layer_3_variant_calling('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/wt_7/', 'wt_7', {'l1': (30, 180), 'l2': (30, 180), 'l3': (30, 180)}, nc=6)

# wt_18
layer_3_variant_calling('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/wt_18/', 'wt_18', {'l1': (40, 200), 'l2': (30, 180), 'l3': (30, 180)}, nc=6)

# wt_19
layer_3_variant_calling('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/wt_19/', 'wt_19', {'l1': (40, 220), 'l2': (40, 220), 'l3': (40, 220)}, nc=6)

# mut_11_1
layer_3_variant_calling('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/mut_11_1/', 'mut_11_1', {'l1': (20, 160), 'l2': (20, 150), 'l3': (20, 180)}, nc=6)

# mut_11_2
layer_3_variant_calling('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/mut_11_2/', 'mut_11_2', {'l1': (30, 200), 'l2': (20, 180), 'l3': (30, 200)}, nc=6)

# mut_15
layer_3_variant_calling('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/mut_15/', 'mut_15', {'l1': (40, 220), 'l2': (30, 220), 'l3': (40, 220)}, nc=6)


def layer_specific_gene_conversion():
    import numpy as np
    import pandas as pd
    from collections import deque
    from tqdm import tqdm
    from matplotlib import pyplot as plt
    import sys
    sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/software/hometools/')
    from myUsefulFunctions import snvdata
    S_COV = {'l1': (20, 150),
             'l2': (20, 130),
             'l3': (20, 180)}
    SAMPLES = ('l1', 'l2', 'l3')
    BASE_DICT = {'A': 4, 'C': 5, 'G': 6, 'T': 7}

    snpsdata = {}
    with open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.snps.txt", 'r') as fin: # Test with lower read depth
        for line in fin:
            line = line.strip().split()
            snpsdata[(line[0], int(line[2]))] = (line[3], line[4])
    # Get read-counts at SNP positions. code in layer_specific_dna_analysis.sh
    # Check Allele frequency distribution at SNP positions
    fig = plt.figure()
    i = 1
    snps = {}
    for sample in SAMPLES:
        af = deque()
        f = f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/mut_11_1_{sample}/{sample}.syri_snps.bamrc'
        with open(f, 'r') as fin:
            for line in tqdm(fin):
                line = line.strip().split()
                q = snpsdata[(line[0], int(line[1]))][1]
                if int(line[3]) == 0:
                    af.append(0)
                    continue
                af.append(round(int(line[BASE_DICT[q]])/int(line[3]), 4))
        ax = fig.add_subplot(3, 1, i)
        ax.hist(af, bins=[v/100 for v in range(101)], range=[0.01, 1])
        ax.set_title(sample)
        ax.grid(which='both', axis='both')
        ax.set_axisbelow(True)
        ax.set_xlabel("Allele Frequency")
        snps[sample] = af
        i += 1
    plt.tight_layout()

    snps = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.snps.txt', header=None)
    snps.drop_duplicates(subset=[0, 1], inplace=True, ignore_index=True)
    for sample in SAMPLES:
        f = f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/mut_11_1_{sample}/{sample}.syri_snps.bamrc'
        af = deque()
        with open(f, 'r') as fin:
            for line in tqdm(fin):
                line = line.strip().split()
                q = snpsdata[(line[0], int(line[1]))][1]
                if int(line[3]) == 0:
                    af.append((line[0], int(line[1]), 0, 0))
                    continue
                af.append((line[0], int(line[1]), round(int(line[BASE_DICT[q]])/int(line[3]), 4), int(line[3])))
        afdict = pd.DataFrame(af)
        afdict.columns = [0, 1, f'{sample}_af', f'{sample}_rc']
        afdict.drop_duplicates(subset=[0, 1], inplace=True, ignore_index=True)
        snps = snps.merge(afdict, how='left', on=[0, 1])
        snps.reset_index(drop=True, inplace=True)

    snps['type'] = 'transversions'
    snps.loc[(snps[3] == 'A') & (snps[4] == 'G'), 'type'] = 'transitions'
    snps.loc[(snps[3] == 'C') & (snps[4] == 'T'), 'type'] = 'transitions'
    snps.loc[(snps[3] == 'G') & (snps[4] == 'A'), 'type'] = 'transitions'
    snps.loc[(snps[3] == 'T') & (snps[4] == 'C'), 'type'] = 'transitions'
    snps.to_csv('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/layers.syri_snps.allele_frequency.txt', sep='\t', index=False)

    # Remove snps where allele frequency is in 0.35-0.6 for all layers
    afdiff = deque()
    for row in snps.itertuples(index=False):
        afdiff.append(max([row[5], row[7], row[9]]) - min([row[5], row[7], row[9]]))
    snps['afdiff'] = afdiff
    snps['mean_rc'] = (snps['l1_rc']+snps['l2_rc']+snps['l3_rc'])/3

    filter = deque()
    for row in snps.itertuples(index=False):
        if 0.35 <= row[5] <= 0.6:
            if 0.35 <= row[7] <= 0.6:
                if 0.35 <= row[9] <= 0.6:
                    filter.append(False)
                    continue
        if row[5] <= 0.1:
            if row[7] <= 0.1:
                if row[9] <= 0.1:
                    filter.append(False)
                    continue
        if max([row[5], row[7], row[9]]) - min([row[5], row[7], row[9]]) < 0.25:
            filter.append(False)
            continue
        filter.append(True)
    snpsfilt = snps.loc[list(filter)].copy()

    ### Filter positions with low sequencing
    snpsfilt = snpsfilt.loc[(snpsfilt['l1_rc'] >= S_COV['l1'][0]) &
                            (snpsfilt['l2_rc'] >= S_COV['l2'][0]) &
                            (snpsfilt['l3_rc'] >= S_COV['l3'][0])].copy()


    def d(rs):
        '''
        This function calculates distance of a point from a line. Here, using it to calculate distance of gene conv candidates from the 3D-diagonal

        Original answer: https://stackoverflow.com/a/50728570
        '''
        p = np.array([0, 0, 0])
        q = np.array([1, 1, 1])
        x = p-q
        return np.linalg.norm(
            np.outer(np.dot(rs-q, x)/np.dot(x, x), x)+q-rs,
            axis=1)
    rs = np.array([[row[5], row[7], row[9]] for row in snpsfilt.itertuples(index=False)])
    dist_diag = d(rs)

    # Select positions for which reads support different haplotypes. Positions for which atleast 5 reads have mismatched haplotypes were selected. Only high quality syri SNP positions are being used

    def get_paired_gene_conv(sample):
        with open(f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/{sample}/paired_snps_gene_conversions.txt', 'w') as fout:
            # with open(f'paired_snps_gene_conversions.txt', 'w') as fout:
            print(sample)
            for i in range(1, 9):
                snps = deque()
                last = -1
                c = f'CUR{i}G'
                with open(f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/{sample}/snps_{c}.pileup', 'r') as fin:
                    # print(sample,c)
                    for line in tqdm(fin):
                        # print(len(snps))
                        line = line.strip().split()
                        try:
                            while True:
                                if int(line[1]) - snps[0].pos > 1000:
                                    snps.popleft()
                                else:
                                    break
                        except IndexError:
                            pass
                        if not S_COV[sample][0] < int(line[3]) < S_COV[sample][1]:
                            continue
                        try:
                            _ = snpsdata[(line[0], line[1])]
                        except KeyError:
                            continue
                        snpd = snvdata(line[:6])
                        setattr(snpd, "BP", list(map(int, line[6].split(','))))
                        reads = line[7].split(',')
                        setattr(snpd, "ref_reads", set([reads[i] for i in range(len(reads)) if snpd.bases[i] in {'.', ','}]))
                        setattr(snpd, "qry_reads", set([reads[i] for i in range(len(reads)) if snpd.bases[i] in {snpsdata[(line[0], line[1])][1], snpsdata[(line[0], line[1])][1].lower()}]))
                        if len(snps) > 0:
                            for s in snps:
                                mismatch = 0          # Different haplotypes
                                rname = deque()
                                if len(snpd.qry_reads) > 0:
                                    for r in s.ref_reads:
                                        if r in snpd.qry_reads:
                                            mismatch += 1
                                            rname.append(r)
                                match = len(s.ref_reads) - mismatch
                                mismatch2 = 0
                                if len(snpd.ref_reads) > 0:
                                    for r in s.qry_reads:
                                        if r in snpd.ref_reads:
                                            mismatch2 += 1
                                            rname.append(r)
                                match += len(s.qry_reads) - mismatch2
                                if mismatch >= 5:
                                    fout.write(f'{s.chr}\t{s.pos}\t{snpd.chr}\t{snpd.pos}\t{mismatch}\t{match}\t{",".join(rname)}\tref_to_qry')
                                    fout.write("\n")
                                if mismatch2 >= 5:
                                    fout.write(f'{s.chr}\t{s.pos}\t{snpd.chr}\t{snpd.pos}\t{mismatch2}\t{match}\t{",".join(rname)}\tqry_to_ref')
                                    fout.write("\n")
                        snps.append(snpd)

    with Pool(processes=8) as pool:
        pool.map(get_paired_gene_conv, SAMPLES)

    # Plot the distribution of read counts with matching and mismatching haplotypes
    from matplotlib.backends.backend_pdf import PdfPages
    with PdfPages(f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/gene_conversion_match_mismatch.pdf') as pdf:
        for sample in SAMPLES:
            with open(f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/{sample}/paired_snps_gene_conversions.txt', 'r') as fin:
                # with open(f'paired_snps_gene_conversions.txt', 'r') as fin:
                mismatch = deque()
                match = deque()
                for line in fin:
                    line = line.strip().split()
                    mismatch.append(int(line[4]))
                    match.append(int(line[5]))
            fig = plt.figure()
            bins = range(0, np.max(match), 5)
            ax1 = plt.subplot2grid((3, 3), (0, 0), rowspan=1, colspan=2)
            ax1.set_xlim([0, np.max(match)+5])
            # ax1.hist(match, bins=int(np.max(match)/5))
            ax1.hist(match, bins=bins)
            ax1.set_ylabel('Frequency')
            ax2 = plt.subplot2grid((3, 3), (1, 0), rowspan=2, colspan=2)
            ax2.set_xlim([0, np.max(match)+5])
            ax2.set_ylim([0, np.max(match)+5])
            ax2.scatter(match, mismatch, s=1, zorder=1, color='black')
            ax2.plot([0, np.max(match)], [0, np.max(match)], zorder=0, label="x1")
            ax2.plot([0, np.max(match)], [0, 0.5*np.max(match)], zorder=0, label="x0.5")
            ax2.plot([0, np.max(match)], [0, 0.25*np.max(match)], zorder=0, label="x0.25")
            ax2.plot([0, np.max(match)], [0, 0.1*np.max(match)], zorder=0, label="x0.1")
            ax2.set_xlabel("# Reads with haplotype match")
            ax2.set_ylabel("# Reads with haplotype mismatch")
            ax2.legend()
            ax3 = plt.subplot2grid((3, 3), (1, 2), rowspan=2, colspan=1)
            ax3.set_ylim([0, np.max(match)+5])
            # ax3.hist(mismatch, orientation = 'horizontal', bins=int(np.max(match)/5))
            ax3.hist(mismatch, orientation='horizontal', bins=bins)
            ax3.set_xlabel('Frequency')
            fig.suptitle(sample, fontsize=12)
            plt.tight_layout()
            pdf.savefig()
            plt.close()

    # Get overlap of selected positions and filter positions that are shared by all samples from a sequencing technology
    genecovs = {sample: pd.read_table(f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/{sample}/paired_snps_gene_conversions.txt', header=None) for sample in SAMPLES}
    genecovpos = {}
    for sample in SAMPLES:
        df = genecovs[sample]
        pos = list(df[0].astype('str') + ':' + df[1].astype('str')) + list(df[2].astype('str') + ':' + df[3].astype('str'))
        genecovpos[sample] = set(pos)

    a = genecovpos['WT_1'].intersection(genecovpos['WT_19']).intersection(genecovpos['MUT_11_1']).intersection(genecovpos['MUT_15'])
    b = genecovpos['wt7'].intersection(genecovpos['wt18']).intersection(genecovpos['mut4']).intersection(genecovpos['mut11_2'])
    badpos = a.union(b)

    genecov = pd.DataFrame()
    for sample in SAMPLES:
        df = genecovs[sample].copy()
        df['sample'] = sample
        df['apos'] = df[0].astype('str') + ':' + df[1].astype('str')
        df['bpos'] = df[2].astype('str') + ':' + df[3].astype('str')
        genecov = pd.concat([genecov, df])
    genecov.reset_index(inplace=True, drop=True)

    goodpos = genecov.loc[~((genecov['apos'].isin(badpos)) | (genecov['bpos'].isin(badpos)))]
    with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/gene_conversions_positions.txt', 'w') as fout:
        allpos = sorted(set(list(goodpos['apos']) + list(goodpos['bpos'])))
        for p in allpos:
            p = p.split(":")
            fout.write(f"{p[0]}:{p[1]}-{p[1]}\n")

    # Get pileup data at these positions. run code in leaf_dna_analysis.sh

    # For the positions selected above, check if reads supporting the same haplotype switch are present in other samples, if yes then remove the positions

    pos_pair = list(zip(list(goodpos['apos']), list(goodpos['bpos'])))
    pos_mismatch_cnt = defaultdict(dict)
    for sample in SAMPLES:
        snvsdata = {}
        with open(f"/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/{sample}/geneconv.pileup", 'r') as fin:
            for line in fin:
                line = line.strip().split()
                if line[3] == '0':
                    continue
                setattr(snpd, "BP", list(map(int, line[6].split(','))))
                snpd = snvdata(line[:6])
                reads = line[7].split(',')
                setattr(snpd, "ref_reads", set([reads[i] for i in range(len(reads)) if snpd.bases[i] in {'.', ','}]))
                setattr(snpd, "qry_reads", set([reads[i] for i in range(len(reads)) if snpd.bases[i] in {snpsdata[(line[0], line[1])][1], snpsdata[(line[0], line[1])][1].lower()}]))
                snvsdata[f"{line[0]}:{line[1]}"] = snpd
        for p in pos_pair:
            try:
                snpd1 = snvsdata[p[0]]
                snpd2 = snvsdata[p[1]]
            except KeyError as e:
                pos_mismatch_cnt[p][sample] = (0, 0, 0)
                continue
            mismatch = 0          # Different haplotypes
            match = 0
            rname = deque()
            if len(snpd2.qry_reads) > 0:
                for r in snpd1.ref_reads:
                    if r in snpd2.qry_reads:
                        mismatch += 1
                        rname.append(r)
            match = len(snpd1.ref_reads) - mismatch
            mismatch2 = 0
            if len(snpd2.ref_reads) > 0:
                for r in snpd1.qry_reads:
                    if r in snpd2.ref_reads:
                        mismatch2 += 1
                        rname.append(r)
            match += len(snpd1.qry_reads) - mismatch2
            # pos_mismatch_cnt[p][sample] = (match, mismatch, mismatch2, rname)
            pos_mismatch_cnt[p][sample] = (match, mismatch, mismatch2)

    select = deque()
    for row in goodpos.itertuples(index=False):
        s = 0
        for k, v in pos_mismatch_cnt[(row[9], row[10])].items():
            if v[1] > 0 or v[2] > 0:
                s += 1
        if s == 8:
            select.append(False)
        else:
            select.append(True)

    goodpos = goodpos[list(select)]
    goodpos.reset_index(inplace=True, drop=True)


    fig = plt.figure(figsize=[8, 4])
    ax = fig.add_subplot(1, 2, 1)
    ax.set_xlim([0, np.max(goodpos[4])])
    ax.hist(goodpos[4], bins=range(np.max(goodpos[4])+1))
    ax.set_xlabel("# Reads with haplotype switch")
    ax.set_ylabel("# Gene conversions")

    ax = fig.add_subplot(1, 2, 2)
    ax.set_xlim([0, 1])
    ax.hist(goodpos[4]/goodpos[5], bins=[i/10 for i in range(11)])
    ax.set_xlabel("AF Reads with haplotype switch")
    ax.set_ylabel("# Gene conversions")
    plt.tight_layout()

    # Most candidate gene conversions are selected because of noisy read mapping
    # To filter them out:
    #       1) Get a bam file with only the reads supporting gene conversions
    #       2) Get pileup data and see if the reads have other linked mutations

    poslist = set(goodpos['apos']).union(goodpos['bpos'])
    for sample in SAMPLES:
        readnames = deque()
        with open(f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/{sample}/geneconv.pileup', 'r') as fin:
            for line in fin:
                line = line.strip().split()
                if f'{line[0]}:{line[1]}' in poslist:
                    readnames.extend(line[7].split(","))
        readnames = list(set(readnames))
        with open(f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/{sample}/geneconv.reads.txt', 'w') as fout:
            for r in readnames:
                fout.write(r +"\t\n")

    import pickle
    # with open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling//gen_conv_goodpos.pickle", 'wb') as fout:
    #     pickle.dump(goodpos, fout)
    with open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling//gen_conv_goodpos.pickle", 'rb') as fin:
        goodpos = pickle.load(fin)

    # Get bam and pileup for the reads in geneconv.reads.txt. code in leaf_dna_analysis.sh
    # Filter out positions where the reads that support gene conversions also have other mutations that are not present in background reads
    def filter_noisy_reads(sample):
        selected = pd.DataFrame()
        readdata = defaultdict(dict)
        posdata = dict()
        with open(f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/{sample}/geneconv.reads.pileup', 'r') as fin:
            for line in tqdm(fin):
                line = line.strip().split()
                pos = snvdata(line[:6])
                if pos.rc != 0:
                    reads = line[7].split(',')
                    setattr(pos, "ref_reads", set([reads[i] for i in range(len(reads)) if pos.bases[i] in {'.', ','}]))
                    setattr(pos, "qry_reads", set([reads[i] for i in range(len(reads)) if pos.bases[i] not in {'.', ','}]))
                    for i in range(len(reads)):
                        readdata[reads[i]][f'{line[0]}:{line[1]}'] = pos.bases[i]
                posdata[f'{line[0]}:{line[1]}'] = pos
                pos.getindelreads(line[4], line[7])

        df = goodpos.loc[goodpos['sample'] == sample]
        s = deque()
        for row in df.itertuples(index=False):
            # if row[1] == 18152399:
            #     break
            bad = False
            gcreads = set(row[6].split(','))
            gcpos = set([k for r in gcreads for k in readdata[r]])
            for p in gcpos:
                # if p in (f"{row[0]}:{row[1]}", f"{row[2]}:{row[3]}"): continue
                gc_m = 0    # Gene conversion reads with reference allele
                gc_mm = 0   # Gene conversion reads with alternate allele
                ref_m = 0
                ref_mm = 0
                gc_ind = 0
                ref_ind = 0
                if sum([1 for r in gcreads if r in posdata[p].qry_reads]) >= 3:
                    if p in (f"{row[0]}:{row[1]}", f"{row[2]}:{row[3]}"): continue
                    #     print(p)
                    for r in posdata[p].ref_reads:
                        if r in gcreads:
                            gc_m += 1
                        else:
                            ref_m += 1
                    for r in posdata[p].qry_reads:
                        if r in gcreads:
                            gc_mm += 1
                        else:
                            ref_mm += 1
                    if gc_mm > ref_mm:
                        bad = True
                        break
                if sum([1 for r in gcreads if r in posdata[p].indelreads]) >= 3:
                    for r in posdata[p].indelreads:
                        if r in gcreads:
                            gc_ind += 1
                        else:
                            ref_ind += 1
                    if gc_ind > ref_ind:
                        bad = True
                        break
            if not bad:
                s.append(row)
        return pd.DataFrame(s)

    with Pool(processes=2) as pool:
        selected = pool.map(filter_noisy_reads, SAMPLES)
    selected = pd.concat(selected)

    # This filtered gene conversion candidates significantly. The remaining candidates still seem to consist of positions that are a result of mis-mapping
    # Sample, candidate GC, filtered GC
    # WT_1 48 8
    # wt7 32 5
    # wt18 20 3
    # WT_19 46 5
    # mut4 29 4
    # MUT_11_1 54 5
    # mut11_2 30 4
    # MUT_15 37 5
    ## When using all SNP_pos
    # WT_1 777 202
    # wt7 316 42
    # wt18 235 34
    # WT_19 1103 267
    # mut4 224 33
    # MUT_11_1 1164 339
    # mut11_2 241 30
    # MUT_15 685 169
    # Next idea is to align the reads to the ORA reference genome and check if they still show same haplotype switch. Assumption is that the most of the mis-mapping reads would be gone when aligning with the other haplotype, whereas true reads supporting true gene-coversion would stay.
    # Code to get reads and align them in leaf_dna_analysis.sh

    # Read the corresponding positions of CUR SNP in the ORA assembly
    selected_snp_pos = set(list(selected.apos) + list(selected.bpos))
    snp_cur_ora_map = {}
    with open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out", 'r') as fin:
        for line in tqdm(fin):
            line = line.strip().split()
            if line[10] != 'SNP': continue
            if f"{line[0]}:{line[1]}" in selected_snp_pos:
                snp_cur_ora_map[f"{line[0]}:{line[1]}"] = f"{line[5]}:{line[6]}"

    sample_concordance = {}
    for sample in SAMPLES:
        readdata = defaultdict(dict)
        posdata = dict()
        with open(f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/{sample}/geneconv.reads.ora.pileup', 'r') as fin:
            for line in fin:
                line = line.strip().split()
                pos = snvdata(line[:6])
                if pos.rc != 0:
                    reads = line[7].split(',')
                    setattr(pos, "ref_reads", set([reads[i] for i in range(len(reads)) if pos.bases[i] in {'.', ','}]))
                    setattr(pos, "qry_reads", set([reads[i] for i in range(len(reads)) if pos.bases[i] not in {'.', ','}]))
                    for i in range(len(reads)):
                        readdata[reads[i]][f'{line[0]}:{line[1]}'] = pos.bases[i]
                posdata[f'{line[0]}:{line[1]}'] = pos
                pos.getindelreads(line[4], line[7])

        concordance = deque()
        for row in selected.loc[selected['sample'] == sample].itertuples(index=False):
            ora1 = snp_cur_ora_map[row.apos]
            ora2 = snp_cur_ora_map[row.bpos]
            c = deque()
            if row[7] == 'ref_to_qry':
                for r in row[6].split(","):
                    try:
                        if readdata[r][ora1] not in [',', '.'] and readdata[r][ora2] in [',', '.']:
                            c.append(True)
                        else:
                            c.append(False)
                    except KeyError:
                        c.append(False)
            elif row[7] == 'qry_to_ref':
                for r in row[6].split(","):
                    try:
                        if readdata[r][ora1] in [',', '.'] and readdata[r][ora2] not in [',', '.']:
                            c.append(True)
                        else:
                            c.append(False)
                    except KeyError:
                        c.append(False)
            concordance.append(sum(c))
        sample_concordance[sample] = concordance

    for sample in SAMPLES:
        print(sample, Counter(sample_concordance[sample]))

#END
layer_specific_gene_conversion()


def test_allele_frequency_variation(cwd, bname, leafrc, step):
    '''
    Layer 3 sequencing has abnormal coverage distribution, suggesting dna degradation.
    Here, I compare the allele frequency at SNP positions in the MUT_11 leaf, layer 1, layer 2, and layer 3.
    If SNPs in L3 have similar allele frequency as other samples, then it would mean that the absurdities in
    coverage are not affecting allele frequencies and that this data can be used for mutation calling.
    '''
    print(cwd, bname, leafrc, step)
    # return
    from tqdm import tqdm
    from collections import deque
    from matplotlib import pyplot as plt
    import numpy as np
    import sys
    import pickle
    plt.interactive(False)
    sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/software/hometools/')
    from myUsefulFunctions import density_scatter, mylogger
    plt.interactive(False)
    logger = mylogger('afplot')

    if step == 1:
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
        leaf_af = {}
        with open(leafrc, 'r') as fin:
            for line in tqdm(fin):
                line = line.strip().split()
                if (line[0], int(line[1])) in snp_pos:
                    leaf_af[(line[0], int(line[1]))] = (int(line[3]), int(line[BD[syri_snp_list[(line[0], int(line[1]))][1]]]))
        layer_af = deque()
        for l in ("l1", "l2", "l3"):
            ldict = {}
            with open(f"{cwd}/{bname}_{l}/bam_read_counts_b30_q10.bt2.txt", "r") as fin:
                for line in tqdm(fin):
                    line = line.strip().split()
                    if (line[0], int(line[1])) in snp_pos:
                        ldict[(line[0], int(line[1]))] = (int(line[3]), int(line[BD[syri_snp_list[(line[0], int(line[1]))][1]]]))
            layer_af.append(ldict)
        with open(f"{cwd}/pickled_snp_af", "wb") as f:
            pickle.dump([leaf_af, layer_af], f)
        return

    elif step == 2:
        with open(f"{cwd}/pickled_snp_af", "rb") as f:
            garb = pickle.load(f)
            leaf_af = garb[0]
            layer_af = garb[1]
        snps_in_all = set(leaf_af.keys()).intersection(set(layer_af[0].keys())).intersection(set(layer_af[1].keys())).intersection(set(layer_af[2].keys()))
        fig = plt.figure(figsize=[12, 12])
        fig.suptitle(bname)
        # plot Leaf vs L3
        ax1 = fig.add_subplot(3, 3, 1)
        ax1.set_xlim([-0.02, 1.02])
        ax1.set_ylim([-0.02, 1.02])
        x = np.array([0 if leaf_af[p][0] == 0 else leaf_af[p][1]/leaf_af[p][0] for p in snps_in_all])
        y = np.array([0 if layer_af[2][p][0] == 0 else layer_af[2][p][1]/layer_af[2][p][0] for p in snps_in_all])
        print(round(np.corrcoef(x, y)[0, 1], 4))
        ax1 = density_scatter(x, y, bins=[[i/100 for i in range(0, 104, 4)], [i/100 for i in range(0, 104, 4)]], ax=ax1, fig=fig, s=0.5)
        ax1.set_xlabel("leaf")
        ax1.set_ylabel("L3")
        ax1.set_title(f"Correlation: {round(np.corrcoef(x, y)[0, 1], 4)}")

        # plot Leaf vs L2
        ax2 = fig.add_subplot(3, 3, 4)
        ax2.set_xlim([-0.02, 1.02])
        ax2.set_ylim([-0.02, 1.02])
        x = np.array([0 if leaf_af[p][0] == 0 else leaf_af[p][1]/leaf_af[p][0] for p in snps_in_all])
        y = np.array([0 if layer_af[1][p][0] == 0 else layer_af[1][p][1]/layer_af[1][p][0] for p in snps_in_all])
        print(round(np.corrcoef(x, y)[0, 1], 4))
        ax2 = density_scatter(x, y, bins=[[i/100 for i in range(0, 104, 4)], [i/100 for i in range(0, 104, 4)]], ax=ax2, fig=fig, s=0.5)
        ax2.set_xlabel("leaf")
        ax2.set_ylabel("L2")
        ax2.set_title(f"Correlation: {round(np.corrcoef(x, y)[0,1], 4)}")

        # plot Leaf vs L1
        ax3 = fig.add_subplot(3, 3, 7)
        ax3.set_xlim([-0.02, 1.02])
        ax3.set_ylim([-0.02, 1.02])
        x = np.array([0 if leaf_af[p][0] == 0 else leaf_af[p][1]/leaf_af[p][0] for p in snps_in_all])
        y = np.array([0 if layer_af[0][p][0] == 0 else layer_af[0][p][1]/layer_af[0][p][0] for p in snps_in_all])
        print(round(np.corrcoef(x, y)[0, 1], 4))
        ax3 = density_scatter(x, y, bins=[[i/100 for i in range(0, 104, 4)], [i/100 for i in range(0, 104, 4)]], ax=ax3, fig=fig, s=0.5)
        ax3.set_xlabel("leaf")
        ax3.set_ylabel("L1")
        ax3.set_title(f"Correlation: {round(np.corrcoef(x, y)[0,1], 4)}")

        # plot L1 vs L3
        ax4 = fig.add_subplot(3, 3, 2)
        ax4.set_xlim([-0.02, 1.02])
        ax4.set_ylim([-0.02, 1.02])
        x = np.array([0 if layer_af[0][p][0] == 0 else layer_af[0][p][1]/layer_af[0][p][0] for p in snps_in_all])
        y = np.array([0 if layer_af[2][p][0] == 0 else layer_af[2][p][1]/layer_af[2][p][0] for p in snps_in_all])
        print(round(np.corrcoef(x, y)[0, 1], 4))
        ax4 = density_scatter(x, y, bins=[[i/100 for i in range(0, 104, 4)], [i/100 for i in range(0, 104, 4)]], ax=ax4, fig=fig, s=0.5)
        ax4.set_xlabel("L1")
        ax4.set_ylabel("L3")
        ax4.set_title(f"Correlation: {round(np.corrcoef(x, y)[0,1], 4)}")

        # plot L1 vs L2
        ax5 = fig.add_subplot(3, 3, 5)
        ax5.set_xlim([-0.02, 1.02])
        ax5.set_ylim([-0.02, 1.02])
        x = np.array([0 if layer_af[0][p][0] == 0 else layer_af[0][p][1]/layer_af[0][p][0] for p in snps_in_all])
        y = np.array([0 if layer_af[1][p][0] == 0 else layer_af[1][p][1]/layer_af[1][p][0] for p in snps_in_all])
        print(round(np.corrcoef(x, y)[0, 1], 4))
        ax5 = density_scatter(x, y, bins=[[i/100 for i in range(0, 104, 4)], [i/100 for i in range(0, 104, 4)]], ax=ax5, fig=fig, s=0.5)
        ax5.set_xlabel("L1")
        ax5.set_ylabel("L2")
        ax5.set_title(f"Correlation: {round(np.corrcoef(x, y)[0,1], 4)}")

        # plot L2 vs L3
        ax6 = fig.add_subplot(3, 3, 3)
        ax6.set_xlim([-0.02, 1.02])
        ax6.set_ylim([-0.02, 1.02])
        x = np.array([0 if layer_af[1][p][0] == 0 else layer_af[1][p][1]/layer_af[1][p][0] for p in snps_in_all])
        y = np.array([0 if layer_af[2][p][0] == 0 else layer_af[2][p][1]/layer_af[2][p][0] for p in snps_in_all])
        print(round(np.corrcoef(x, y)[0, 1], 4))
        ax6 = density_scatter(x, y, bins=[[i/100 for i in range(0, 104, 4)], [i/100 for i in range(0, 104, 4)]], ax=ax6, fig=fig, s=0.5)
        ax6.set_xlabel("L2")
        ax6.set_ylabel("L3")
        ax6.set_title(f"Correlation: {round(np.corrcoef(x, y)[0, 1], 4)}")
        plt.tight_layout()
        plt.savefig(f"{cwd}/SNP_AF_correlation.png", dpi=300)

    else:
        logger.error("Incorrect value for step. Required values: 1 or 2.")
# END
import os
from multiprocessing import Pool
from functools import partial
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
samples = ('wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_15')
cwds = [cwd + 'layer_samples/' + s + os.sep for s in samples]
leafrcs = ['WT_1/bam_read_counts_b30_q10.bt2.txt',
           'wt7/bam_read_counts_b30_q10.bt2.txt',
           'wt18/bam_read_counts_b30_q10.bt2.txt',
           'WT_19/bam_read_counts_b30_q10.bt2.txt',
           'MUT_11_1/bam_read_counts_b30_q10.bt2.txt',
           'mut11_2/bam_read_counts_b30_q10.bt2.txt',
           'MUT_15/bam_read_counts_b30_q10.bt2.txt']
leafrcs = [cwd+s for s in leafrcs]
params = list(zip(cwds, samples, leafrcs))
with Pool(processes=8) as pool:
    pool.starmap(partial(test_allele_frequency_variation, step=1), params)
with Pool(processes=8) as pool:
    pool.starmap(partial(test_allele_frequency_variation, step=2), params)


############################# Sub-functions ####################################

def getsyrivarlist():
    syriout = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out'
    syri_snp_list = {}
    syri_indel_list = {}
    with open(syriout, 'r') as fin:
        # types = {'SNP', 'INS', 'DEL'}
        for line in fin:
            line = line.strip().split()
            if line[10] == 'SNP':
                if line[3] == 'N' or line[4] == 'N': continue
                syri_snp_list[(line[0], int(line[1]))] = (line[3], line[4])
            elif line[10] == 'INS':
                syri_indel_list[(line[0], int(line[1]))] = ('+', line[4][1:])
            elif line[10] == 'DEL':
                syri_indel_list[(line[0], int(line[1])+1)] = ('-', line[3][1:])
    return syri_snp_list, syri_indel_list
# END


def readfilteredbamreadcount(filename, scov_min, scov_max, noise_pos, syri_snp_list, syri_indel_list):
    from collections import deque
    from tqdm import tqdm

    sample_noise = deque()
    snp_alfrq = {}
    indel_alfrq = {}
    sm_pos = {}
    snp_pos = set(syri_snp_list.keys())
    indel_pos = set(syri_indel_list.keys())
    base_dict = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
    base_dict2 = {4: 'A', 5: 'C', 6: 'G', 7: 'T'}

    with open(filename, 'r') as fin:
        for line in tqdm(fin):
            line = line.strip().split()
            # Get allele frequency if position is in syri_snp
            if (line[0], int(line[1])) in snp_pos:
                p = (line[0], int(line[1]))
                snp_alfrq[p] = (round(int(line[base_dict[syri_snp_list[p][1]]])/int(line[3]), 3), int(line[3]))
                continue
            # Get allele frequency if position is in syri_indel
            if (line[0], int(line[1])) in indel_pos:
                p = (line[0], int(line[1]))
                v = syri_indel_list[p][0] + syri_indel_list[p][1]
                for j in range(9, len(line), 2):
                    if line[j] == v:
                        indel_alfrq[p] = (round(int(line[j+1])/int(line[3]), 3), int(line[3]))
                        break
                continue
            # Check if the position is noisy, if not then select somatic mutation
            if (line[0], int(line[1])) not in noise_pos:
                # Check read-depth at the position
                if not scov_min <= int(line[3]) <= scov_max:
                    sample_noise.append((line[0], int(line[1])))
                    continue
                ind = [4, 5, 6, 7] + list(range(10, len(line), 2))
                # for i in ind:
                rc = 0
                ar = 0
                alt = ''
                for i in ind:
                    if int(line[i]) >= 5:
                        if base_dict[line[2]] != i:
                            if int(line[i]) <= rc:
                                continue
                            rc = int(line[i])
                            ar = round(int(line[i])/int(line[3]), 3)
                            try:
                                alt = base_dict2[i]
                            except KeyError:
                                alt = line[i-1]
                if rc >= 5:
                    sm_pos[(line[0], int(line[1]))] = (line[2], alt, rc, ar)
    return sample_noise, snp_alfrq, indel_alfrq, sm_pos
# END
