import os
from collections import defaultdict, OrderedDict, deque
cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
SAMPLES = ("MUT_11_1", "MUT_15", "WT_1", "WT_19")
BASE_DICT = {'A': 4,
             'C': 5,
             'G': 6,
             'T': 7}

for sample in SAMPLES:
    candidates = {}
    goodcand = {}
    with open('{cwd}/{sample}/{sample}_good_candidates.txt'.format(cwd=cwd, sample=sample)) as fin:
        for line in fin:
            line = line.strip().split()
            candidates[line[0] + '_' + line[2]] = (line[3], line[4], line[5], line[6])

    with open('{cwd}/{sample}/{sample}_good_candidates_rna_read_counts.txt'.format(cwd=cwd, sample=sample)) as fin:
        for line in fin:
            line = line.strip().split()
            if int(line[BASE_DICT[candidates[line[0] + '_' + line[1]][1]]]) >= 2:
                print(line)
    break







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