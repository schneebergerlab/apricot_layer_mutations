from collections import deque, Counter, defaultdict
from tqdm import tqdm
import numpy as np
import pandas as pd
from copy import deepcopy
import sys
sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/software/hometools/')
from myUsefulFunctions import snvdata
from matplotlib import pyplot as plt
from multiprocessing import Pool

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
indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling'
SAMPLES = ('WT_1', 'wt7', 'wt18', 'WT_19', 'mut4', 'MUT_11_1', 'mut11_2', 'MUT_15')
BASE_DICT = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
BASE_DICT2 = {4: 'A', 5: 'C', 6: 'G', 7: 'T'}
S_COV = {'WT_1': (50, 230),
         'wt7': (40, 150),
         'wt18': (40, 150),
         'WT_19': (50, 240),
         'mut4': (40, 140),
         'MUT_11_1': (50, 250),
         'mut11_2': (40, 140),
         'MUT_15': (50, 230)}
# Testing with lower read depth
S_COV = {'WT_1': (20, 230),
         'wt7': (20, 150),
         'wt18': (20, 150),
         'WT_19': (20, 240),
         'mut4': (20, 140),
         'MUT_11_1': (20, 250),
         'mut11_2': (20, 140),
         'MUT_15': (20, 230)}
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
    with open(f'{indir}/{sample}/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt') as fin:
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

## Remove positions that are present in all samples
## Will do this based on sequencing technologies
to_pop = deque()
for k in samples_sm['WT_1']:
    try:
        for sample in ('WT_19', 'MUT_11_1', 'MUT_15'):
            _ = samples_sm[sample][k]
        to_pop.append(k)
    except KeyError:
        pass

for k in samples_sm['wt7']:
    try:
        for sample in ('wt18', 'mut4', 'mut11_2'):
            _ = samples_sm[sample][k]
        to_pop.append(k)
    except KeyError:
        pass

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
    ax = fig.add_subplot(4, 4, i)
    ax.hist([v[3] for v in samples_sm[sample].values()], range=[0, 1], bins=100, color=SAMCOLS[sample])
    ax.set_yscale('log')
    ax.set_title(f'{sample} allele frequency')
    ax.axvline(x=0.25, color='black', ls=':')
    i += 1
    ax = fig.add_subplot(4, 4, i)
    ax.hist([v[2] for v in samples_sm[sample].values()], range=[0, 100], bins=100, color=SAMCOLS[sample])
    ax.set_yscale('log')
    ax.set_title(f'{sample} read count')
    ax.axvline(x=20, color='black', ls=':')
    i += 1
plt.tight_layout()
plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/af_rc_after_initial_filtering.pdf')
plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/af_rc_after_initial_filtering.png')

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
plot_selected_pos(hcpos, "tmp_igv.bat", "tmp/hcsm/", M=75, HEIGHT=600, emptydir=True)


# Write hcsm to table and then check individual positions manually to select high confidence mutations
CWD = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
posd = deque()
for k, l in hcsm.items():
    for pos, v in l.items():
        posd.append(list(pos) + list(v) + [k])
posd = pd.DataFrame(posd)
posd.sort_values([0, 1], inplace=True)
posd.to_csv(f"{CWD}/high_cov_mutants_sorted.all_samples.txt", sep=',', index=False, header=False)


# Manually curated list of somatic mutations in the 8 branches
high_cov_mutants_sorted.all_samples.txt

# Remove hcsm from sample_sm and then use the resulting dict for filtering low conf somatic mutations
for sample in SAMPLES:
    for k in hcsm[sample]:
        samples_sm[sample].pop(k)

# First filter low conf SM using the filtered positions allele frequency and then using the raw-data from bam-read-count output file

### Remove noisy positions that are in the samples_sm object
to_pop = deque()
for sample in SAMPLES:
    for k, v in tqdm(samples_sm[sample].items()):
        # Filter out positions that have >0.2x reads in the background samples.
        rcnt = 0
        for s2 in SAMPLES:
            try:
                v_s2 = samples_sm[s2][k]
            except KeyError:
                continue
            rcnt += v_s2[2]
            if rcnt > 1.2*v[2]:
                to_pop.append(k)

for sample in SAMPLES:
    print(sample, len(samples_sm[sample]))
    for p in to_pop:
        try:
            samples_sm[sample].pop(p)
        except KeyError:
            pass
    print(sample, len(samples_sm[sample]))

### Remove noisy positions
SMP_SM_BLP = deepcopy(samples_sm)
def remove_noisy_positions(fin, SAMPLES, samples_sm):
    tpop = set()
    pos = set([k for s in SAMPLES for k in samples_sm[s]])
    pos_sam = {k: s for s in SAMPLES for k in samples_sm[s]}
    # pos_rcnt = {p: 0 for p in pos}
    for s in SAMPLES:
        spop = deque()
        with open(f'{indir}/{s}/{fin}') as f:
            for line in tqdm(f):
                line = line.strip().split()
                p = (line[0], int(line[1]))
                if p not in pos:
                    continue
                if p in tpop:
                    continue
                if s == pos_sam[p]:
                    continue
                ## Filter positions that have > 1 or more than 0.2x alt allele read count in other samples compared to focal sample
                # scnt = 0.2*samples_sm[pos_sam[p]][p][2]
                # try:
                #     ind = BASE_DICT[samples_sm[pos_sam[p]][p][1]]
                #     if int(line[ind]) > 1:
                #         spop.append(p)
                #         continue
                #     if int(line[ind]) >= scnt:
                #         spop.append(p)
                #         continue
                # except KeyError:
                #     v = samples_sm[pos_sam[p]][p][1]
                #     for i in range(9, len(line), 2):
                #         if line[i] == v:
                #             if int(line[i+1]) > 1:
                #                 spop.append(p)
                #                 continue
                #             if int(line[i+1]) >= scnt:
                #                 spop.append(p)
                #                 continue
                ind = [4, 5, 6, 7] + list(range(10, len(line), 2))
                ind.remove(BASE_DICT[line[2]])
                scnt = int(0.2*samples_sm[pos_sam[p]][p][2])
                rcnt = 0
                for i in ind:
                    # pos_rcnt[(line[0], p)] += int(line[i])
                    rcnt += int(line[i])
                if rcnt >= scnt:
                    spop.append(p)
        tpop = tpop.union(set(spop))
    return tpop

# from filtered read count
to_pop = remove_noisy_positions('filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', SAMPLES, samples_sm)
for sample in SAMPLES:
    for p in to_pop:
        try:
            samples_sm[sample].pop(p)
        except KeyError:
            pass

# from bam read count
# TODO: filter using b0_q0 read counts
## tried to parallelize this, but multiprocessing would not work because each process tries to copy all objects in memory, resulting in too much memory usage. Multithreading uses less memory but is not faster than normal.
to_pop = remove_noisy_positions('bam_read_counts_b0_q0.bt2.txt', SAMPLES, samples_sm)
for sample in SAMPLES:
    for p in to_pop:
        try:
            samples_sm[sample].pop(p)
        except KeyError:
            pass


## Gene conversion identification at SNP positions
## Plot allele frequency at the heterozygous positions

fig = plt.figure(figsize=[10, 6])
i = 1
snps = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.snps.txt', header=None)

for sample in SAMPLES:
    f = f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/{sample}/{sample}.syri_snps.bamrc'
    af = deque()
    with open(f, 'r') as fin:
        for line in tqdm(fin):
            line = line.strip().split()
            q = syri_snp_list[(line[0], int(line[1]))][1]
            if int(line[3]) == 0:
                af.append(0)
                continue
            af.append(round(int(line[BASE_DICT[q]])/int(line[3]), 4))
    ax = fig.add_subplot(4, 2, i)
    ax.hist(af, bins=100, range=[0.01, 1])
    ax.set_title(sample)
    ax.grid(which='both', axis='both')
    ax.set_axisbelow(True)
    snps[sample] = af
    i += 1
plt.tight_layout()
plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_samples.syri_snps.allele_frequency.pdf')
plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_samples.syri_snps.allele_frequency.png')

snps = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.snps.txt', header=None)
snps.drop_duplicates(subset=[0, 1], inplace=True, ignore_index=True)
for sample in SAMPLES:
    f = f'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/{sample}/{sample}.syri_snps.bamrc'
    af = deque()
    with open(f, 'r') as fin:
        for line in tqdm(fin):
            line = line.strip().split()
            q = syri_snp_list[(line[0], int(line[1]))][1]
            if int(line[3]) == 0:
                af.append((line[0], int(line[1]), 0, 0))
                continue
            af.append((line[0], int(line[1]), round(int(line[BASE_DICT[q]])/int(line[3]), 4), int(line[3])))
    afdict = pd.DataFrame(af)
    afdict.columns = [0, 1, f'{sample}_af', f'{sample}_rc']
    afdict.drop_duplicates(subset=[0, 1], inplace=True, ignore_index=True)
    snps = snps.merge(afdict, how='left', on=[0, 1])
    snps.reset_index(drop=True, inplace=True)
snps.to_csv('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_samples.syri_snps.allele_frequency.txt', sep='\t', index=False)

### Filter positions with low sequencing
snpsfilt = snps.loc[(snps['wt7_rc'] >= 20) &
                    (snps['wt18_rc'] >= 20) &
                    (snps['mut4_rc'] >= 20) &
                    (snps['mut11_2_rc'] >= 20)].copy()

snpsfilt = snpsfilt.loc[(snpsfilt['wt7_af'] <= 0.7) &
                        (snpsfilt['wt18_af'] <= 0.7) &
                        (snpsfilt['WT_1_af'] <= 0.7) &
                        (snpsfilt['WT_19_af'] <= 0.7) &
                        (snpsfilt['wt7_af'] >= 0.25) &
                        (snpsfilt['wt18_af'] >= 0.25) &
                        (snpsfilt['WT_1_af'] >= 0.25) &
                        (snpsfilt['WT_19_af'] >= 0.25)].copy()

snpsfilt = snpsfilt.loc[((snpsfilt['mut4_af'] >= 0.95) |
                        (snpsfilt['mut11_2_af'] >= 0.95) |
                        (snpsfilt['MUT_11_1_af'] >= 0.95) |
                        (snpsfilt['MUT_15_af'] >= 0.95)) |
                        ((snpsfilt['mut4_af'] <= 0.05) |
                        (snpsfilt['mut11_2_af'] <= 0.05) |
                        (snpsfilt['MUT_11_1_af'] <= 0.05) |
                        (snpsfilt['MUT_15_af'] <= 0.05))]

from matplotlib.backends.backend_pdf import PdfPages
with PdfPages("gene_conversion_candidates_af.pdf") as pdf:
    count = 0
    try:
        del fig
    except:
        pass
    for row in snpsfilt.itertuples(index=False):
        if count % 10 == 0:
            try:
                plt.tight_layout()
                pdf.savefig(fig)
                plt.close()
                # break
            except NameError:
                pass
            count = 0
            fig = plt.figure(figsize=[10, 14])
        count += 1
        ax = fig.add_subplot(10, 1, count)
        ax.set_ylim([0, 1])
        ax.set_xlim([0, 7])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.plot(list(range(8)), [row[i] for i in range(5, 21, 2)], linewidth=2)
        ax.set_xticklabels(SAMPLES)
        ax.set_title(f'{row[0]}:{row[1]}-{row[2]}')
        ax.minorticks_on()
        ax.grid(which='both', axis='y')
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()
## Cannot find any good gene conversion candidates

## Checking only the new samples
# Select positions for which average AF in WT  is >0.3 and <0.7, where average AF in MUT in <0.15 and >0.85
gene_cov = dict()
for k, v in tqdm(syri_snp_list.items()):
    good = False
    wt_af = 0
    mut_af = 0
    # for sample in ['WT_1', 'WT_19']:
    for sample in ['wt7', 'wt18']:
        try:
            wt_af = sample_syri_snp[sample][k][0]
            wt_rc = sample_syri_snp[sample][k][1]
        except KeyError:
            # continue
            good = False
            break
        if 0.2 <= wt_af <= 0.8:
            if S_COV[sample][0] <= wt_rc <= S_COV[sample][1]:
                good = True
                # break
            else:
                # continue
                good = False
                break
        else:
            # continue
            good = False
            break
    if not good:
        continue
    # for sample in ['MUT_11_1', 'MUT_15']:
    good = False
    for sample in ['mut4',  'mut11_2']:
        try:
            mut_af = sample_syri_snp[sample][k][0]
            mut_rc = sample_syri_snp[sample][k][1]
        except KeyError:
            good = True
            continue
            # gene_cov[k] = v
            # good = False
            # break
        if not S_COV[sample][0] <= mut_rc <= S_COV[sample][1]:
            good = False
            break
            # continue
        if mut_af > 0.05 and mut_af < 0.95:
            good = False
            # good = True
            break
    if good:
        gene_cov[k] = v

k = list(gene_cov.keys())
a = np.random.randint(0, len(gene_cov), 20)
gcpos = set([f'{i[0]}_{i[1]}' for i in sorted(gene_cov.keys())])
gcpos = set([f'{i[0]}_{i[1]}' for i in sorted([k[i] for i in a])])
plot_selected_pos(gcpos, "tmp_igv.bat", "tmp/gene_conv/", M=75, HEIGHT=600, emptydir=True)

#-------------------------------------------------------------------------------
## Gene conversion identification by comparing/using linked SNP positions
snps = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.snps.txt', header=None)
snps[1] = snps[1].astype(int)
i = 1
# Get distribution of inter-snp distance
fig = plt.figure()
for grp in snps.groupby(0):
    garb = grp[1].sort_values([1])
    x = np.array(garb.loc[:, 1])[1:]
    y = np.array(garb.loc[:, 1])[:-1]
    print(grp[0], min(x-y), max(x-y))
    ax = fig.add_subplot(4, 2, i)
    ax.hist(x-y, bins=100, range=[0, 1000])
    # ax.hist([grp[1][1][i] - grp[1][1][i-1] for i in range(1, grp[1].shape[0])])
    ax.set_title(grp[0])
    ax.set_yscale('log')
    i += 1

# Get all snps that are within 1kb to each other
snpsfilt = pd.DataFrame()
for grp in snps.groupby(0):
    garb = grp[1].sort_values([1])
    # x = np.array(garb[1][1:])
    # y = np.array(garb[1][:-1])
    x = np.array(garb.loc[:, 1])[1:]
    y = np.array(garb.loc[:, 1])[:-1]
    diff = x-y < 1000
    close_snps = garb.loc[[diff[0]] + list(diff)]
    snpsfilt = pd.concat([snpsfilt, close_snps])
snpsfilt[1] -= 1
snpsfilt.to_csv('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/snps_close.bed', index=False, sep='\t', header=None)

### Get linked indel positions as well (adding for layer_enriched data, leaf data can also be analysed later)
snps = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.snps.txt', header=None)
snps.columns = ['chromosome', 'start', 'end', 'ref_allele', 'alt_allele']
indels = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.indels.txt', header=None)
indels.columns = ['chromosome', 'start', 'end', 'alt_allele']
shv = pd.concat([snps, indels])
shv['start'] = shv['start'].astype(int)
shv.sort_values(['chromosome', 'start'], inplace=True)

i = 1
# Get distribution of inter-snp distance
fig = plt.figure()
for grp in shv.groupby('chromosome'):
    garb = grp[1].sort_values(['start'])
    x = np.array(garb.loc[:, 'start'])[1:]
    y = np.array(garb.loc[:, 'start'])[:-1]
    print(grp[0], min(x-y), max(x-y))
    ax = fig.add_subplot(4, 2, i)
    ax.hist(x-y, bins=100, range=[0, 1000])
    # ax.hist([grp[1][1][i] - grp[1][1][i-1] for i in range(1, grp[1].shape[0])])
    ax.set_title(grp[0])
    ax.set_yscale('log')
    i += 1

# Get all snps that are within 1kb to each other
shvsfilt = pd.DataFrame()
for grp in shv.groupby('chromosome'):
    garb = grp[1].sort_values(['start'])
    # x = np.array(garb[1][1:])
    # y = np.array(garb[1][:-1])
    x = np.array(garb.loc[:, 'start'])[1:]
    y = np.array(garb.loc[:, 'start'])[:-1]
    diff = x-y < 1000
    close_shv = garb.loc[[diff[0]] + list(diff)]
    shvsfilt = pd.concat([shvsfilt, close_shv])
shvsfilt['start'] -= 1
shvsfilt.to_csv('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/shv_close.bed', index=False, sep='\t', header=None)


# Get pileup data at SNP positions (done in leaf_dna_analysis.sh)
import sys
sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/software/hometools/')
from myUsefulFunctions import snvdata
snpsdata = {}
# with open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/strict_syn_snp_allele_readcount.depth350-550.af0.35-0.6.bed", 'r') as fin:
with open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/snps_close.bed", 'r') as fin: # Test with lower read depth
    for line in fin:
        line = line.strip().split()
        snpsdata[(line[0], line[2])] = (line[3], line[4])


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

"""
After running the above test, 3 positions have shown consistent signal for gene conversion (1 in WT_1, and 2 in WT_19). However, these are 10x libraries with PCR amplified DNA and as such these gene conversions are doubtful.
Sample:Concordant gene conversion reads in CUR and ORA:Number of reads supporting gene conversion (as identified when using CUR as ref)  
WT_1 : deque([1, 0, 0, 0, 0, 0, 5, 2]) [5, 9, 5, 7, 5, 5, 5, 8]
wt7 : deque([0, 0, 0, 0, 0]) [5, 7, 7, 7, 7]
wt18 : deque([0, 0, 0]) [5, 5, 5]
WT_19 : deque([1, 5, 6, 0, 0]) [6, 5, 6, 5, 6]
mut4 : deque([0, 0, 0, 0]) [6, 6, 6, 8]
MUT_11_1 : deque([0, 1, 0, 0, 0]) [5, 5, 5, 5, 6]
mut11_2 : deque([0, 0, 0, 0]) [6, 5, 5, 5]
MUT_15 : deque([0, 0, 0, 0, 0]) [5, 6, 6, 5, 5]

# With all SNPs
WT_1 Counter({0: 170, 1: 15, 2: 6, 6: 4, 5: 3, 3: 2, 7: 1, 4: 1})
wt7 Counter({0: 42})
wt18 Counter({0: 28, 1: 6})
WT_19 Counter({0: 228, 1: 13, 5: 11, 3: 5, 2: 5, 4: 3, 7: 1, 6: 1})
mut4 Counter({0: 32, 5: 1})
MUT_11_1 Counter({0: 283, 1: 26, 5: 10, 4: 7, 2: 4, 6: 4, 3: 3, 7: 1, 12: 1})
mut11_2 Counter({0: 29, 1: 1})
MUT_15 Counter({0: 128, 4: 10, 1: 9, 5: 8, 2: 5, 6: 4, 7: 3, 3: 2}) 
"""
