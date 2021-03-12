import os
from collections import defaultdict, OrderedDict, deque
cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
SAMPLES = ("MUT_11_1", "MUT_15", "WT_1", "WT_19")
BASE_DICT = {'A': 0,
             'C': 1,
             'G': 2,
             'T': 3
             }

snps_readcount = dict()

i = 0
for sample in SAMPLES:
    print(sample)
    sample_dict = defaultdict()
    os.chdir('{}/{}'.format(cwd, sample))
    with open('filtered_low_ref_al_bam_read_counts_b30_q10.txt'.format(sample), 'r') as fin:
        for line in fin:
            line = line.strip().split()
            sample_dict[line[0] + '_' + line[1]] = (line[4], line[5], line[6], line[7])
    snps_readcount[sample] = sample_dict

candidate_readcount = dict()
for sample in SAMPLES:
    print(sample)
    sample_dict = defaultdict()
    os.chdir('{}/{}'.format(cwd, sample))
    with open('multi_cell_{}_filtered_SNPs_candidate.sorted.bed'.format(sample), 'r') as fin:
        for line in fin:
            line = line.strip().split()
            sample_dict[line[0] + '_' + line[2]] = (line[3], line[4], int(line[5]))
    candidate_readcount[sample] = sample_dict

candidate_readcount = OrderedDict(candidate_readcount)
import matplotlib.pyplot as plt

pcnt = 1
fig = plt.figure(figsize=[10,10])
for sample in SAMPLES:
    sample_readcount = defaultdict(deque)
    for pos, v in candidate_readcount[sample].items():
        for s in SAMPLES:
            try:
                sample_readcount[s].append(int(snps_readcount[s][pos][BASE_DICT[v[1]]]))
            except KeyError as e:
                sample_readcount[s].append(0)
        # print(sample, chr, pos, v[2])
    # print(sample_readcount)
    ax = fig.add_subplot(2,2,pcnt)
    pcnt += 1
    marker = ['+', 'x', '.']
    m = 0
    for s in SAMPLES:
        if s == sample: continue
        ax.scatter(sample_readcount[sample], sample_readcount[s], label=s, marker=marker[m])
        m+=1
    ax.legend()
    ax.set_xlabel('{} Alt AF'.format(sample))
    ax.set_ylabel('Other sample Alt AF')
    ax.set_ylim([-0.5, 30])
    ax.set_xlim([-0.5, 30])
plt.tight_layout()
# plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/alt_allele_readcounts.png')
plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/alt_allele_readcounts_zoomed.png')
plt.close()