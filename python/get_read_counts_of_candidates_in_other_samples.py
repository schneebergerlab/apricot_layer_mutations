import os
from collections import defaultdict, OrderedDict, deque
cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
SAMPLES = ("MUT_11_1", "MUT_15", "WT_1", "WT_19")
BASE_DICT = {'A': 0,
             'C': 1,
             'G': 2,
             'T': 3
             }

candidate_readcount = dict()
for sample in SAMPLES:
    print(sample)
    sample_dict = defaultdict()
    os.chdir('{}/{}'.format(cwd, sample))
    with open('multi_cell_bt2_{}_bt2_filtered_SNPs_candidate.sorted.bed'.format(sample), 'r') as fin:
        for line in fin:
            line = line.strip().split()
            sample_dict[line[0] + '_' + line[2]] = (line[3], line[4], int(line[5]), float(line[6]))
    candidate_readcount[sample] = sample_dict
candidate_readcount = OrderedDict(candidate_readcount)

loc = set([k for sample in SAMPLES for k in candidate_readcount[sample].keys()])

from multiprocessing import Pool
from functools import partial

def get_readcounts(fin, loc):
    sample_dict = {}
    with open(fin, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            if line[0] + '_' + line[1] in loc:
                sample_dict[line[0] + '_' + line[1]] = (int(line[4]), int(line[5]), int(line[6]), int(line[7]))
    return sample_dict

snps_readcount = []
fins = []
for sample in SAMPLES:
    fin='{}{}/bam_read_counts_b30_q10.bt2.txt'.format(cwd, sample)
    # fin='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/{}/filtered_low_ref_al_bam_read_counts_b30_q10.txt'.format(sample)
    fins.append(fin)

with Pool(processes = 4) as pool:
    snps_readcount = pool.map(partial(get_readcounts, loc=loc), fins)
import matplotlib.pyplot as plt
snps_readcount = dict(zip(SAMPLES, snps_readcount))

for sample in SAMPLES:
    # Get alt readcount at candidate positions in all samples
    sample_readcount = defaultdict(deque)
    for pos, v in candidate_readcount[sample].items():
        for s in SAMPLES:
            try:
                sample_readcount[s].append(int(snps_readcount[s][pos][BASE_DICT[v[1]]]))
            except KeyError as e:
                # sample_readcount[s].append(0)
                print('Readcount at position: {} not found in sample: {}'.format(pos, sample))     # Use when using the entire bam_readcount file

    # Filter positions where the 'other' samples have as many alt reads as the focal samples
    bad_pos = deque()
    for i in range(len(sample_readcount[sample])):
        s_rc = sample_readcount[sample][i]
        s_all_rc = sum([sample_readcount[s][i] for s in SAMPLES])
        if s_all_rc >= 2*s_rc:
            bad_pos.append(i)
    bad_pos = set(bad_pos)

    # Select positions where the number of alt reads in other samples is less than
    # or equal to 20% of the alt reads count in focal sample
    good_pos = list(set(range(len(sample_readcount[sample]))) - bad_pos)
    high_conf = deque()
    for i in good_pos:
        s_rc = sample_readcount[sample][i]
        s_all_rc = sum([sample_readcount[s][i] for s in SAMPLES])
        if s_all_rc <= 1.2*s_rc:
            high_conf.append(i)


    # Get alt allele-freq at candidate positions in all samples
    sample_alfreq = defaultdict(deque)
    for pos, v in candidate_readcount[sample].items():
        for s in SAMPLES:
            try:
                sample_alfreq[s].append((snps_readcount[s][pos][BASE_DICT[v[1]]])/sum(snps_readcount[s][pos]))
            except KeyError as e:
                # sample_readcount[s].append(0)
                print('Readcount at position: {} not found in sample: {}'.format(pos, sample))     # Use when using the entire bam_readcount file
            except ZeroDivisionError as e:
                sample_alfreq[s].append(1)

    # Filter positions where the cumulative alt allele freq is more than or equal
    # to the focal sample
    bad_pos_af = deque()
    for i in range(len(sample_alfreq[sample])):
        s_rc = sample_alfreq[sample][i]
        s_all_rc = sum([sample_alfreq[s][i] for s in SAMPLES])
        if s_all_rc >= 2*s_rc:
            bad_pos_af.append(i)
    bad_pos_af = set(bad_pos_af)

    # Select positions where the alt allele-freq in other samples is less than or
    # equal to the alt allele freq in the focal sample
    good_pos_af = list(set(range(len(sample_alfreq[sample]))) - bad_pos_af)
    high_conf_af = deque()
    for i in good_pos_af:
        s_rc = sample_alfreq[sample][i]
        s_all_rc = (sum([sample_alfreq[s][i] for s in SAMPLES]) - s_rc)/3
        if s_all_rc <= 0.2*s_rc:
            high_conf_af.append(i)

    all_bad  = set(list(bad_pos) + list(bad_pos_af))
    all_high = set([i for i in high_conf if i in high_conf_af])
    all_mid  = set(range(len(sample_alfreq[sample]))) - set(list(all_bad) + list(all_high))

    with open('{}{}/{}_uncertain_candidates.txt'.format(cwd, sample, sample), 'w') as fout:
        keys = list(candidate_readcount[sample].keys())
        for i in sorted(list(all_mid)):
            chr, p = keys[i].rsplit('_', 1)
            d = [str(c) for c in candidate_readcount[sample][keys[i]]]
            fout.write('\t'.join([chr, str(int(p)-1), p] + d) + '\n')

    with open('{}{}/{}_good_candidates.txt'.format(cwd, sample, sample), 'w') as fout:
        keys = list(candidate_readcount[sample].keys())
        for i in sorted(list(all_high)):
            chr, p = keys[i].rsplit('_', 1)
            d = [str(c) for c in candidate_readcount[sample][keys[i]]]
            fout.write('\t'.join([chr, str(int(p)-1), p] + d) + '\n')

    positions_set = [set(range(len(sample_readcount[sample]))), all_bad, all_mid, all_high]
    pcnt = 1
    fig = plt.figure(figsize=[10, 10])
    marker = ['+', 'x', '.']
    plt.rc('font', size = 8)
    for i in range(4):
        positions = positions_set[i]
        ax = fig.add_subplot(4, 2, pcnt)
        pcnt += 1
        m = 0
        for s in SAMPLES:
            if s == sample: continue
            ax.scatter([sample_readcount[sample][i] for i in positions], [sample_readcount[s][i] for i in positions], label=s, marker=marker[m])
            m+=1
        ax.legend(loc='upper left')
        ax.set_xlabel('{}'.format(sample))
        ax.set_ylabel('Other samples')
        if i == 0: title, lim = ('Alt allele readcount: All Candidates (#{})'.format(len(positions)), 250)
        elif i == 1: title, lim = ('Alt allele readcount: Bad Candidates (#{})'.format(len(positions)), 250)
        elif i == 2: title, lim = ('Alt allele readcount: Uncertain Candidates (#{})'.format(len(positions)), 250)
        elif i == 3: title, lim = ('Alt allele readcount: Good Candidates (#{})'.format(len(positions)), 250)
        ax.set_title(title)
        ax.set_ylim([-1, lim + 1])
        ax.set_xlim([-1, lim + 1])

        ax = fig.add_subplot(4, 2, pcnt)
        pcnt += 1
        m = 0
        for s in SAMPLES:
            if s == sample: continue
            ax.scatter([sample_alfreq[sample][i] for i in positions], [sample_alfreq[s][i] for i in positions], label=s, marker=marker[m])
            m+=1
        ax.legend(loc='upper left')
        ax.set_xlabel('{}'.format(sample))
        ax.set_ylabel('Other samples')
        if i == 0: title, lim = ('Alt allele frequency: All Candidates (#{})'.format(len(positions)), 1.05)
        elif i == 1: title, lim = ('Alt allele frequency: Bad Candidates (#{})'.format(len(positions)), 1.05)
        elif i == 2: title, lim = ('Alt allele frequency: Uncertain Candidates (#{})'.format(len(positions)), 1.05)
        elif i == 3: title, lim = ('Alt allele frequency: Good Candidates (#{})'.format(len(positions)), 1.05)
        ax.set_title(title)
        ax.set_ylim([-0.05, lim])
        ax.set_xlim([-0.05, lim])
    plt.tight_layout()
    plt.savefig('{}{}/{}_alt_alleles.png'.format(cwd, sample, sample))
    plt.close()

