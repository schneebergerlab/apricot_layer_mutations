import os
from collections import defaultdict, OrderedDict, deque
from tqdm import tqdm
from multiprocessing import Pool
from functools import partial
import pandas as pd
from subprocess import run
import sys
sys.path.insert(0, '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python')
from get_candidate_for_mutant_phenotype_funcs import getlinefrombamrc, tepid_te_del_igv_batch, tefinder_igv_batch
from multiprocessing import set_start_method, Pool
# set_start_method('spawn')
from matplotlib import pyplot as plt
import numpy as np
from datetime import datetime
from upsetplot import UpSet, from_contents
import pybedtools as bdt
from pybedtools.cbedtools import Interval
import pyranges as pr

CWD = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/'
INDIR = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
SAMPLES = ("MUT_11_1", "MUT_15", "WT_1", "WT_19")
BASE_DICT = {'A': 0,
             'C': 1,
             'G': 2,
             'T': 3
             }

############## This is the logical pipeline, but the results are weird #########
# SNP Identification
candidate_readcount = dict()
for sample in SAMPLES:
    sample_dict = defaultdict()
    os.chdir('{}/{}'.format(INDIR, sample))
    with open('multi_cell_bt2_{}_bt2_filtered_SNPs_candidate.sorted.bed'.format(sample), 'r') as fin:
        for line in fin:
            line = line.strip().split()
            sample_dict[line[0] + '_' + line[2]] = (line[3], line[4], int(line[5]), float(line[6]))
    candidate_readcount[sample] = sample_dict
candidate_readcount = OrderedDict(candidate_readcount)

wt1pos = set(candidate_readcount['WT_1'].keys())
wt19pos = set(candidate_readcount['WT_19'].keys())
mut11pos = set(candidate_readcount['MUT_11_1'].keys())
mut15pos = set(candidate_readcount['MUT_15'].keys())
mut_pos = set([k for k in mut11pos if k in mut15pos])
locmt = set([m for m in mut_pos if candidate_readcount['MUT_11_1'][m][:2] == candidate_readcount['MUT_15'][m][:2]])
wt_pos = set([k for k in wt1pos if k in wt19pos])
locwt = set([m for m in wt_pos if candidate_readcount['WT_1'][m][:2] == candidate_readcount['WT_19'][m][:2]])
loc = set(locwt).union(locmt)

def get_readcounts(fin, loc):
    sample_dict = {}
    with open(fin, 'r') as fin:
        for line in tqdm(fin):
            line = line.strip().split()
            if line[0] + '_' + line[1] in loc:
                sample_dict[line[0] + '_' + line[1]] = (int(line[4]), int(line[5]), int(line[6]), int(line[7]))
    return sample_dict
snps_readcount = []
fins = []
for sample in SAMPLES:
    fin='{}{}/bam_read_counts_b30_q10.bt2.txt'.format(INDIR, sample)
    fins.append(fin)
with Pool(processes=4) as pool:
    snps_readcount = pool.map(partial(get_readcounts, loc=loc), fins)
snps_readcount = dict(zip(SAMPLES, snps_readcount))

# Select positions where the number of alt reads in WT samples is less than
# or equal to 20% of the alt reads count in MUT sample with lower alt read-count
goodmt = deque()
for l in tqdm(sorted(locmt)):
    refi = BASE_DICT[candidate_readcount['MUT_11_1'][l][0]]
    alti = BASE_DICT[candidate_readcount['MUT_11_1'][l][1]]
    altrc = min(candidate_readcount['MUT_11_1'][l][2], candidate_readcount['MUT_15'][l][2])
    if snps_readcount['WT_1'][l][alti] >= 0.2*altrc or snps_readcount['WT_19'][l][alti] >= 0.2*altrc:
        continue
    altfr = min(candidate_readcount['MUT_11_1'][l][3], candidate_readcount['MUT_15'][l][3])
    if sum(snps_readcount['WT_19'][l]) == 0 or sum(snps_readcount['WT_1'][l]) == 0: continue
    if snps_readcount['WT_1'][l][alti]/sum(snps_readcount['WT_1'][l]) >= 0.2*altfr or snps_readcount['WT_19'][l][alti]/sum(snps_readcount['WT_19'][l]) >= 0.2*altfr:
        continue
    goodmt.append(l.rsplit("_", maxsplit=1))

# Filter mutations that have another SNP within 100bp
chrs = set([i[0] for i in goodmt])
chrdic = defaultdict(deque)
for v in goodmt:
    chrdic[v[0]].append(int(v[1]))
for k, v in chrdic.items():
    vs = sorted(v)
    bp = deque()
    for i in range(1, len(vs)):
        if vs[i] - vs[i-1] < 100:
            bp.append(vs[i-1])
            bp.append(vs[i])
    bp = set(bp)
    gp = sorted(list(set(vs) - set(bp)))
    chrdic[k] = gp
goodmt = [k+"_"+str(p) for k, v in chrdic.items() for p in v]


with open(CWD + "denovo_mutation.txt", 'w') as fout:
    fout.write("chr{t}pos{t}ref{t}alt{t}wt1_rc{t}wt1_af{t}wt19_rc{t}wt19_af{t}mut11_rc{t}mut11_af{t}mut15_rc{t}mut15_af".format(t="\t") + "\n")
    for l in goodmt:
        refi = BASE_DICT[candidate_readcount['MUT_11_1'][l][0]]
        alti = BASE_DICT[candidate_readcount['MUT_11_1'][l][1]]
        fout.write("\t".join(l.rsplit("_", 1) +
                             list(candidate_readcount['MUT_11_1'][l][:2]) +
                             [str(snps_readcount['WT_1'][l][alti]), str(round(snps_readcount['WT_1'][l][alti]/sum(snps_readcount['WT_1'][l]), 4))] +
                             [str(snps_readcount['WT_19'][l][alti]), str(round(snps_readcount['WT_19'][l][alti]/sum(snps_readcount['WT_19'][l]), 4))] +
                             [str(snps_readcount['MUT_11_1'][l][alti]), str(round(snps_readcount['MUT_11_1'][l][alti]/sum(snps_readcount['MUT_11_1'][l]), 4))] +
                             [str(snps_readcount['MUT_15'][l][alti]), str(round(snps_readcount['MUT_15'][l][alti]/sum(snps_readcount['MUT_15'][l]), 4))]) + "\n")


goodwt = deque()
for l in tqdm(sorted(locwt)):
    refi = BASE_DICT[candidate_readcount['WT_1'][l][0]]
    alti = BASE_DICT[candidate_readcount['WT_1'][l][1]]
    # Select positions where the number of alt reads in WT samples is less than
    # or equal to 20% of the alt reads count in MUT sample with lower alt read-count
    altrc = min(candidate_readcount['WT_1'][l][2], candidate_readcount['WT_19'][l][2])
    if snps_readcount['MUT_11_1'][l][alti] >= 0.2*altrc or snps_readcount['MUT_15'][l][alti] >= 0.2*altrc:
        continue
    altfr = min(candidate_readcount['WT_1'][l][3], candidate_readcount['WT_19'][l][3])
    if sum(snps_readcount['MUT_11_1'][l]) == 0 or sum(snps_readcount['MUT_15'][l]) == 0: continue
    if snps_readcount['MUT_11_1'][l][alti]/sum(snps_readcount['MUT_11_1'][l]) >= 0.2*altfr or snps_readcount['MUT_15'][l][alti]/sum(snps_readcount['MUT_15'][l]) >= 0.2*altfr:
        continue
    goodwt.append(l.rsplit("_", maxsplit=1))

# Filter mutations that have another SNP within 100bp
chrs = set([i[0] for i in goodwt])
chrdic = defaultdict(deque)
for v in goodwt:
    chrdic[v[0]].append(int(v[1]))
for k,v in chrdic.items():
    vs = sorted(v)
    bp = deque()
    for i in range(1, len(vs)):
        if vs[i] - vs[i-1] < 100:
            bp.append(vs[i-1])
            bp.append(vs[i])
    bp = set(bp)
    gp = sorted(list(set(vs) - set(bp)))
    chrdic[k] = gp
goodwt = [k+"_"+str(p) for k, v in chrdic.items() for p in v]

with open(CWD + "mutation_lose.txt", 'w') as fout:
    fout.write("chr{t}pos{t}ref{t}alt{t}wt1_rc{t}wt1_af{t}wt19_rc{t}wt19_af{t}mut11_rc{t}mut11_af{t}mut15_rc{t}mut15_af".format(t="\t") + "\n")
    for l in goodwt:
        refi = BASE_DICT[candidate_readcount['WT_1'][l][0]]
        alti = BASE_DICT[candidate_readcount['WT_1'][l][1]]
        fout.write("\t".join(l.rsplit("_", 1) +
                             list(candidate_readcount['WT_1'][l][:2]) +
                             [str(snps_readcount['WT_1'][l][alti]), str(round(snps_readcount['WT_1'][l][alti]/sum(snps_readcount['WT_1'][l]), 4))] +
                             [str(snps_readcount['WT_19'][l][alti]), str(round(snps_readcount['WT_19'][l][alti]/sum(snps_readcount['WT_19'][l]), 4))] +
                             [str(snps_readcount['MUT_11_1'][l][alti]), str(round(snps_readcount['MUT_11_1'][l][alti]/sum(snps_readcount['MUT_11_1'][l]), 4))] +
                             [str(snps_readcount['MUT_15'][l][alti]), str(round(snps_readcount['MUT_15'][l][alti]/sum(snps_readcount['MUT_15'][l]), 4))]) + "\n")

# Filter positions for which other samples get >5 reads when no mapping quality filter is used
# OR which have <3 reads when mapping quality > 40

# Merge the target positions to a bam-readcount suitable input
os.chdir(CWD)
pos = deque()
with open("denovo_mutation.txt", 'r') as fin:
    line = fin.readline()
    for line in fin:
        line = line.strip().split()
        pos.append([line[0], line[1], line[1], line[2], line[3], 'MT'])
        # fout.write("\t".join([line[0], line[1], line[1], line]) + "\n")

with open("mutation_lose.txt", 'r') as fin:
    line = fin.readline()
    for line in fin:
        line = line.strip().split()
        pos.append([line[0], line[1], line[1], line[2], line[3], 'WT'])
        # fout.write("\t".join([line[0], line[1], line[1]]) + "\n")
pos = pd.DataFrame(pos)
pos[[1, 2]] = pos[[1, 2]].astype(int)
pos.sort_values([0, 1, 2], inplace=True)
pos.reset_index(inplace=True, drop=True)
pos.to_csv(CWD+"snp_candidates.txt", sep='\t', header=False, index=False)

for sample in SAMPLES:
    os.chdir(CWD)
    # f = open("TMP.txt", 'w')
    run("/srv/netscratch/dep_mercier/grp_schneeberger/bin/hometools pbamrc -w 0 -b 0 -q 0 -l snp_candidates.txt -f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta -n 12 ../{s}/{s}.sorted.bt2.bam {s}_snp_readcount.b0.q0.txt".format(s=sample).split())

    run("/srv/netscratch/dep_mercier/grp_schneeberger/bin/hometools pbamrc -w 0 -b 30 -q 40 -l snp_candidates.txt -f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta -n 12 ../{s}/{s}.sorted.bt2.bam {s}_snp_readcount.b30.q40.txt".format(s=sample).split())

# Filter positions for which other samples get >5 reads when no mapping quality filter is used
# OR which have <3 reads when mapping quality > 40
dfnomap = {}
dfhmap = {}
BASE_DICT = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
for sample in SAMPLES:
    d = pd.read_table("{}/{}_snp_readcount.b0.q0.txt".format(CWD, sample), header=None, engine='python', index_col=False)
    d[9] = sample
    dfnomap[sample] = d

    d = pd.read_table("{}/{}_snp_readcount.b30.q40.txt".format(CWD, sample), header=None, engine='python', index_col=False)
    d[9] = sample
    dfhmap[sample] = d

isbad = deque()
sgrp = {'MT': ('MUT_11_1', 'MUT_15'), 'WT': ('WT_1', 'WT_19')}
for i in range(pos.shape[0]):
    s = pos.iat[i, 5]
    if dfhmap[sgrp[s][0]].iat[i, BASE_DICT[pos.iat[i, 4]]] < 3:
        isbad.append(True)
        continue
    if dfhmap[sgrp[s][1]].iat[i, BASE_DICT[pos.iat[i, 4]]] < 3:
        isbad.append(True)
        continue
    ac = 0
    for sample in SAMPLES:
        if sample in sgrp[s]: continue
        ac += dfnomap[sample].iat[i, BASE_DICT[pos.iat[i, 4]]]
    if ac >= 5:
        isbad.append(True)
        continue
    isbad.append(False)

gooddf = pos.iloc[[i for i in range(len(isbad)) if isbad[i] is False]].copy()
# gooddf.to_csv(CWD+'all_good_candidate.qualfilter.txt', index=False, header=False, sep=' ')
gooddf[1] -= 1
gooddf.to_csv(CWD+'snp_candidates.qualfilter.bed', index=False, header=False, sep=' ')
gooddf[1] += 1

def samplempileup(sample, CWD):
    os.chdir(CWD)
    run("/srv/netscratch/dep_mercier/grp_schneeberger/software/TMP/samtools-1.13/bin/samtools view -@ 3 -b -o {s}_snp_candidates.qualfilter.bam -L snp_candidates.qualfilter.bed ../{s}/{s}.sorted.bt2.bam".format(s=sample).split())
    with open('{s}_snp_candidates.qualfilter.baq.bam'.format(s=sample), 'wb') as f:
        run("/srv/netscratch/dep_mercier/grp_schneeberger/software/TMP/samtools-1.13/bin/samtools calmd -Abr {s}_snp_candidates.qualfilter.bam /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta".format(s=sample).split(), stdout=f)
    run("/srv/netscratch/dep_mercier/grp_schneeberger/software/TMP/samtools-1.13/bin/samtools index {s}_snp_candidates.qualfilter.baq.bam".format(s=sample).split())
    try:
        os.remove('{s}_snp_candidates.qualfilter.baq.mpileup'.format(s=sample))
    except FileNotFoundError:
        pass
    with open('{s}_snp_candidates.qualfilter.baq.mpileup'.format(s=sample), 'w') as fout:
        for row in gooddf.itertuples(index=False):
            run("/srv/netscratch/dep_mercier/grp_schneeberger/software/TMP/samtools-1.13/bin/samtools mpileup -aa -A -f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta -q 40 -Q 13 -r {}:{}-{} {s}_snp_candidates.qualfilter.baq.bam".format(row[0], row[1], row[2], s=sample).split(), stdout=fout)

with Pool(processes=4) as pool:
    out = pool.map(partial(samplempileup, CWD=CWD), SAMPLES)


def getbases(l):
    indelcount = 0
    bases = deque()
    skip = 0
    indel = False
    for c in l:
        if skip > 0 and indel == False:
            skip -= 1
            continue
        if indel == True:
            if c in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                skip = (skip*10) + int(c)
                continue
            # skip = int(c)
            else:
                indel = False
                skip -= 1
                continue
        if c == '*':
            # self.indelcount += 1
            indelcount += 1
            bases.append(c)
            continue
        if c == '$': continue
        if c in ['<', '>']: print('spliced alignment found')
        if c == '^':
            skip = 1
            continue
        if c in ['+', '-']:
            indel = True
            continue
        bases.append(c)
    return [indelcount, list(bases)]


isgood = [True]*gooddf.shape[0]
for sample in SAMPLES:
    os.chdir(CWD)
    mpileup = pd.read_table('{s}_snp_candidates.qualfilter.baq.mpileup'.format(s=sample), header=None, index_col=False)
    for i in range(gooddf.shape[0]):
        if sample in sgrp[gooddf.iat[i, 5]]:
            _, bases = getbases(mpileup.iat[i, 4])
            if len([c for c in bases if c in [gooddf.iat[i, 4].upper(), gooddf.iat[i, 4].lower()]]) < 5:
                isgood[i] = False


low_cov_mutants = gooddf.iloc[isgood].copy()
low_cov_mutants.to_csv(CWD+"/low_cov_mutants.txt", header=False, index=False, sep='\t')
# The selected 3 candidates were manually curated to get high confidence low coverage somatic mutations

################################################################################
################# Test after merging the reads from two mutants ################
################################################################################
CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/'
mut_cand = {}
with open(CWD+'mutant.filtered.regions', 'r') as fin:
    for line in fin:
        line = line.strip().split()
        mut_cand[line[0]+"_"+line[2]] = [(line[3], line[4]), int(line[5]), round(float(line[6]), 4)]
ks = set(list(mut_cand.keys()))
posinwt1 = deque()
with open("../WT_1/WT_1_bt2_SNPs_candidate.sorted.bed", 'r') as fin:
    for line in fin:
        line = line.strip().split()
        if line[0] + "_" + line[2] in ks:
            if (line[3], line[4]) == mut_cand[line[0] + "_" + line[2]][0]:
                posinwt1.append(line[0] + "_" + line[2])

posinwt19 = deque()
with open("../WT_19/WT_19_bt2_SNPs_candidate.sorted.bed", 'r') as fin:
    for line in fin:
        line = line.strip().split()
        if line[0] + "_" + line[2] in ks:
            if (line[3], line[4]) == mut_cand[line[0] + "_" + line[2]][0]:
                posinwt19.append(line[0] + "_" + line[2])

BASE_DICT = {'A': 4,
             'C': 5,
             'G': 6,
             'T': 7}

candks = ks - set(list(posinwt1) + list(posinwt19))
with open("candks.bed", 'w') as fout:
    for i in candks:
        c, p = i.rsplit("_", 1)
        fout.write("\t".join([c, p, p]) + "\n")

def runbamrc(c):
    print(c)
    run(c.split())
cmd = "hometools pbamrc -n 10 -b 30 -q 10 -w 0 -f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta -l candks.bed ../{s}/{s}.sorted.bt2.bam {s}.candks.rc.txt"
cmds = [cmd.format(s=i) for i in SAMPLES]
with Pool(processes=4) as pool:
    pool.map(runbamrc, cmds)



fins = ["../{}/bam_read_counts_b30_q10.bt2.txt".format(i) for i in SAMPLES]
with Pool(processes=4) as pool:
    bamdf = pool.map(partial(getlinefrombamrc, cand=candks), fins)
bamdfmt11 = bamdf[0]
bamdfmt15 = bamdf[1]
bamdfwt1 = bamdf[2]
bamdfwt19 = bamdf[3]

noisewt1 = deque()
for line in bamdfwt1:
    if int(line[BASE_DICT[mut_cand[line[0] + "_" + line[1]][0][1]]])/int(line[3]) >= 0.2*mut_cand[line[0] + "_" + line[1]][2]:
        noisewt1.append(line[0] + "_" + line[1])

noisewt19 = deque()
for line in bamdfwt19:
    if int(line[BASE_DICT[mut_cand[line[0] + "_" + line[1]][0][1]]])/int(line[3]) >= 0.2*mut_cand[line[0] + "_" + line[1]][2]:
        noisewt19.append(line[0] + "_" + line[1])

candks2 = candks - set(noisewt1).union(set(noisewt19))
##
garb = set([k for k in candks if mut_cand[k][2] >= 0.1 and mut_cand[k][1] >= 10])
garbdf = {k: {"MUT": mut_cand[k][1]} for k in garb}
for k, v in {"WT_1": bamdfwt1, "WT_19": bamdfwt19, "MUT_11": bamdfmt11, "MUT_15": bamdfmt15}.items():
    for p in v:
        if p[0] + '_' + p[1] not in garb:
            continue
        pos = p[0] + '_' + p[1]
        garbdf[pos][k] = int(p[BASE_DICT[mut_cand[pos][0][1]]])
garbdf = pd.DataFrame(garbdf).transpose()
garbdf['id'] = garbdf.index.values
garbdf['mt11_ratio'] = 1 - abs((2*(garbdf['MUT_11']/garbdf['MUT'])) - 1)

import plotly.express as px
fig = px.scatter_3d(garbdf, z='MUT', x='WT_1', y='WT_19', color='mt11_ratio', opacity=0.7, size_max=2, hover_data=['WT_1', 'WT_19', 'MUT_11', 'MUT_15', 'MUT'])
fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
fig.show(renderer="browser")

import pickle
with open("garbdf.pickle", 'wb') as fout:
    pickle.dump(garbdf, fout)
with open("garbdf.pickle", 'rb') as fout:
    garbdf = pickle.load(fout)


# Test candidates without any quality-filters
CWD = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/'
BASE_DICT = {'A': 4,
             'C': 5,
             'G': 6,
             'T': 7}
BASE_DICT2 = {'A': 0,
             'C': 1,
             'G': 2,
             'T': 3}
BASES = 'ACGT'
mut_cand = {}
with open(CWD+'mutant.b0.q0.filtered.txt', 'r') as fin:
    for line in tqdm(fin):
        line = line.strip().split()
        if int(line[3]) < 100: continue
        bprc = (int(line[4]), int(line[5]), int(line[6]),  int(line[7]))
        h = [i for i in [0, 1, 2, 3] if bprc[i] >= 50]
        if len(h) != 2: continue
        for i in h:
            if i == BASE_DICT2[line[2]]: continue
            if bprc[i] >= 50:
                altbp = BASES[i]
        mut_cand[line[0]+"_"+line[1]] = [line[2], altbp, bprc]
ks = set(list(mut_cand.keys()))
highks = ks

fins = ["{}.b0.q0.txt".format(i) for i in SAMPLES]
with Pool(processes=4) as pool:
    bamdf = pool.map(partial(getlinefrombamrc, cand=highks), fins)

bamdfmt11 = bamdf[0]
bamdfmt15 = bamdf[1]
bamdfwt1 = bamdf[2]
bamdfwt19 = bamdf[3]

noisewt1 = deque()
for p in bamdfwt1:
    if p[0] + '_' + p[1] in highks:
        pos = p[0] + '_' + p[1]
        if int(p[BASE_DICT[mut_cand[pos][1]]]) >= 15:
            noisewt1.append(pos)

noisewt19 = deque()
for p in bamdfwt19:
    if p[0] + '_' + p[1] in highks:
        pos = p[0] + '_' + p[1]
        if int(p[BASE_DICT[mut_cand[pos][1]]]) >= 15:
            noisewt19.append(pos)
noisewt = set(noisewt1).intersection(set(noisewt19))

candks = highks - noisewt
garbdf = {k: {"MUT": mut_cand[k][2][BASE_DICT2[mut_cand[k][1]]]} for k in candks}
for k, v in {"WT_1": bamdfwt1, "WT_19": bamdfwt19, "MUT_11": bamdfmt11, "MUT_15": bamdfmt15}.items():
    for p in tqdm(v):
        if p[0] + '_' + p[1] not in candks:
            continue
        pos = p[0] + '_' + p[1]
        garbdf[pos][k] = int(p[BASE_DICT[mut_cand[pos][1]]])
garbdf = pd.DataFrame(garbdf).transpose()
garbdf['id'] = garbdf.index.values
garbdf['mt11_ratio'] = [round(i, 3) for i in (1 - abs((2*(garbdf['MUT_11']/garbdf['MUT'])) - 1))]

garbdf2 = garbdf.loc[(garbdf['MUT'] <= 1000)]
import plotly.express as px

fig = px.scatter_3d(garbdf2, z='MUT', x='WT_1', y='WT_19', color='mt11_ratio', opacity=0.7, size_max=2, hover_data=['WT_1', 'WT_19', 'MUT_11', 'MUT_15', 'MUT'])
fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
fig.show(renderer="browser")

# Allele Frequency
garbdf = {k: {"MUT": round(mut_cand[k][2][BASE_DICT2[mut_cand[k][1]]]/sum(mut_cand[k][2]), 3)} for k in candks}
for k, v in {"WT_1": bamdfwt1, "WT_19": bamdfwt19, "MUT_11": bamdfmt11, "MUT_15": bamdfmt15}.items():
    for p in tqdm(v):
        if p[0] + '_' + p[1] not in candks:
            continue
        pos = p[0] + '_' + p[1]
        garbdf[pos][k] = int(p[BASE_DICT[mut_cand[pos][1]]])/sum([int(p[4]),int(p[5]),int(p[6]), int(p[7])])
garbdf = pd.DataFrame(garbdf).transpose()
garbdf['id'] = garbdf.index.values
garbdf['mt11_ratio'] = [round(i, 3) for i in (1 - abs((2*(garbdf['MUT_11']/garbdf['MUT'])) - 1))]
garbdf['mt11_ratio'] = [round(i, 3) for i in (1 - abs((2*(garbdf['MUT_11']/garbdf['MUT'])) - 1))]

# garbdf2 = garbdf.loc[(garbdf['MUT'] <= 1000) & (garbdf['WT_1'] <= 25) & (garbdf['WT_19'] <= 25)]
garbdf2 = garbdf.loc[(garbdf['MUT'] <= 1000)]
import plotly.express as px

fig = px.scatter_3d(garbdf2, z='MUT', x='WT_1', y='WT_19', color='mt11_ratio', opacity=0.7, size_max=2, hover_data=['WT_1', 'WT_19', 'MUT_11', 'MUT_15', 'MUT'])
fig.update_layout(margin=dict(l=0, r=0, b=0, t=0))
fig.show(renderer="browser")

## Compare the unfiltered bam readcount for the four samples directly
w1 = open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/WT_1.b0.q0.txt", 'r')
w2 = open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/WT_19.b0.q0.txt", 'r')
m1 = open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/MUT_11_1.b0.q0.txt", 'r')
m2 = open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/MUT_15.b0.q0.txt", 'r')
BASE_DICT = {'A': 0,
             'C': 1,
             'G': 2,
             'T': 3
             }
posdata = defaultdict(dict)
READSIZE = 100000
goodpos = {}
files = {'w1': w1, 'w2': w2, 'm1': m1, 'm2': m2}
curchr = {'w1': deque(), 'w2': deque(), 'm1': deque(), 'm2': deque()}
cnt = 0
while True:
    print(cnt, str(datetime.now()))
    cnt+=1
    # if cnt == 5: break
    changed = False
    # Get allele counts at positions for all samples
    for k, v in files.items():
        for i in range(READSIZE):
            if w1l == '': break
            else:
                w1l = v.readline().strip().split()
                changed = True
                x = list(map(int, w1l[4:8]))
                posdata[(w1l[0], w1l[1], w1l[2])][k] = x
                if len(curchr[k]) == 0: curchr[k].append(w1l[0])
                elif curchr[k][-1] != w1l[0]: curchr[k].append(w1l[0])

    rem = deque()
    for k, v in posdata.items():
        if len(v) < 4:
            missc = [c for c in curchr.keys() if c not in list(v.keys())]
            for mc in missc:
                if curchr[mc][-1] != k[0]:
                    rem.append(k)
                    break
            continue
        rem.append(k)
        if k[-1] == 'N': continue
        skip = BASE_DICT[k[-1]]

        # Remove candidates where multiple alt alleles have >= 5 reads or where no sample has >=5 alt reads
        hashi = False
        for s, c in v.items():
            high = sum([True if c1 >= 5 else False for c1 in c])
            if high > 2: break
            if high == 2: hashi = True
        if high > 2: continue
        if not hashi: continue

        # Remove candidates where samples have different alleles with >= 5 reads
        m = []
        for s, c in v.items():
            for i in range(4):
                if c[i] >= 5:
                    if i == skip: continue
                    m.append(i)
                    break
        if len(set(m)) > 1: continue
        # if len(set(m)) == 0: print("ERROR in finding mutation position")
        m = m[0]

        # Get minor allele count and frequency in the four samples, and select cand which have 2 fold difference between WT and MUT
        wtrc = v['w1'][m] + v['w2'][m]
        wtaf = wtrc/sum(v['w1'] + v['w2'])
        mtrc = v['m1'][m] + v['m2'][m]
        mtaf = mtrc/sum(v['m1'] + v['m2'])
        if wtaf < 0.5*mtaf: goodpos[k] = v
        elif wtaf > 2*mtaf: goodpos[k] = v
    for r in rem: posdata.pop(r)
    if not changed: break

import pickle
# with open(CWD+"goodpos.pickle", 'wb') as fout:
#     pickle.dump(goodpos, fout)
with open(CWD+"goodpos.pickle", 'rb') as fin:
    goodpos = pickle.load(fin)
pos = list(goodpos.keys())
BASE_DICT = {'A': 0,
             'C': 1,
             'G': 2,
             'T': 3
             }
posdata = {}
for p in tqdm(pos):
    m = []
    skip = BASE_DICT[p[2]]
    for s, c in goodpos[p].items():
        for i in range(4):
            if c[i] >= 5:
                if i == skip: continue
                m.append(i)
                break
        if len(m) > 0:
            m = m[0]
            break
    # Checking only those positions that have at least alt allele frequency of 0.1 in the WT or the MUT branches
    if (goodpos[p]['w1'][m] + goodpos[p]['w2'][m])/(sum(goodpos[p]['w1']) + sum(goodpos[p]['w2'])) > 0.1:
        posdata[p] = (goodpos[p]['w1'][m] + goodpos[p]['w2'][m],
                      sum(goodpos[p]['w1']) + sum(goodpos[p]['w2']),
                      goodpos[p]['m1'][m] + goodpos[p]['m2'][m],
                      sum(goodpos[p]['m1']) + sum(goodpos[p]['m2']))
        continue
    if (goodpos[p]['m1'][m] + goodpos[p]['m2'][m])/(sum(goodpos[p]['m1']) + sum(goodpos[p]['m2'])) > 0.1:
        posdata[p] = (goodpos[p]['w1'][m] + goodpos[p]['w2'][m],
                      sum(goodpos[p]['w1']) + sum(goodpos[p]['w2']),
                      goodpos[p]['m1'][m] + goodpos[p]['m2'][m],
                      sum(goodpos[p]['m1']) + sum(goodpos[p]['m2']))
        continue

fig = plt.figure(figsize=[6,6])
ax = fig.add_subplot(2,1,1)
ax.scatter([v[0] for v in posdata.values()], [v[2] for v in posdata.values()], s=0.5)
ax.set_xlim(-10, 200)
ax.set_ylim(-10, 200)
ax.set_title("Alternate allele read-count")
ax.set_xlabel("Merged WT")
ax.set_ylabel("Merged MUT")

ax = fig.add_subplot(2,1,2)
ax.scatter([v[0]/v[1] for v in posdata.values()], [v[2]/v[3] for v in posdata.values()], s=0.5)
ax.set_xlim(-0.05, 1)
ax.set_ylim(-0.05, 1)
ax.set_title("Alternate allele frequency")
ax.set_xlabel("Merged WT")
ax.set_ylabel("Merged MUT")
plt.tight_layout()
plt.savefig(CWD+"mut_vs_wt_alt_all_rc_af.no_qual_filt.pdf")
plt.savefig(CWD+"mut_vs_wt_alt_all_rc_af.no_qual_filt.png", transparent=True)

## Check crosslabelling of WT and MUT samples
candidate_readcount = dict()
for sample in SAMPLES:
    sample_dict = defaultdict()
    os.chdir('{}/{}'.format(INDIR, sample))
    with open('multi_cell_bt2_{}_bt2_filtered_SNPs_candidate.sorted.bed'.format(sample), 'r') as fin:
        for line in fin:
            line = line.strip().split()
            sample_dict[line[0] + '_' + line[2]] = (line[3], line[4], int(line[5]), float(line[6]))
    candidate_readcount[sample] = sample_dict
candidate_readcount = OrderedDict(candidate_readcount)

### Check if any of the WT is a MUT
wt1pos = set(candidate_readcount['WT_1'].keys())
wt19pos = set(candidate_readcount['WT_19'].keys())
mut11pos = set(candidate_readcount['MUT_11_1'].keys())
mut15pos = set(candidate_readcount['MUT_15'].keys())
mut_pos = mut11pos.intersection(mut15pos)
locmt = set([m for m in mut_pos if candidate_readcount['MUT_11_1'][m][:2] == candidate_readcount['MUT_15'][m][:2]])
locwt19 = locmt.intersection(wt19pos) - wt1pos
locwt19 = set([m for m in locwt19 if candidate_readcount['MUT_11_1'][m][:2] == candidate_readcount['WT_19'][m][:2]])
locwt1 = locmt.intersection(wt1pos) - wt19pos
locwt1 = set([m for m in locwt1 if candidate_readcount['MUT_11_1'][m][:2] == candidate_readcount['WT_1'][m][:2]])
loc = set(locwt1).union(locwt19)
def get_readcounts(fin, loc):
    sample_dict = {}
    with open(fin, 'r') as fin:
        for line in tqdm(fin):
            line = line.strip().split()
            if line[0] + '_' + line[1] in loc:
                sample_dict[line[0] + '_' + line[1]] = (int(line[4]), int(line[5]), int(line[6]), int(line[7]))
    return sample_dict
snps_readcount = []
fins = []
for sample in SAMPLES:
    fin='{}{}/bam_read_counts_b30_q10.bt2.txt'.format(INDIR, sample)
    fins.append(fin)
with Pool(processes=4) as pool:
    snps_readcount = pool.map(partial(get_readcounts, loc=loc), fins)
snps_readcount = dict(zip(SAMPLES, snps_readcount))
goodw19 = deque()
for s in locwt19:
    refi = BASE_DICT[candidate_readcount['WT_19'][s][0]]
    alti = BASE_DICT[candidate_readcount['MUT_11_1'][s][1]]
    altrc = np.mean([candidate_readcount[p][s][2] for p in ['MUT_11_1', 'MUT_15', 'WT_19']])
    altaf = sum([candidate_readcount[p][s][2] for p in ['MUT_11_1', 'MUT_15', 'WT_19']])/sum([sum(snps_readcount[p][s]) for p in ['MUT_11_1', 'MUT_15', 'WT_19']])
    bgrc = snps_readcount['WT_1'][s][alti]
    bgaf = snps_readcount['WT_1'][s][alti]/sum(snps_readcount['WT_1'][s])
    if bgrc <= 0.5*altrc and bgaf <= 0.5*altaf:
        goodw19.append((s, altrc, altaf, bgrc, bgaf))
goodw1 = deque()
for s in locwt1:
    refi = BASE_DICT[candidate_readcount['WT_1'][s][0]]
    alti = BASE_DICT[candidate_readcount['WT_1'][s][1]]
    altrc = np.mean([candidate_readcount[p][s][2] for p in ['MUT_11_1', 'MUT_15', 'WT_1']])
    altaf = sum([candidate_readcount[p][s][2] for p in ['MUT_11_1', 'MUT_15', 'WT_1']])/sum([sum(snps_readcount[p][s]) for p in ['MUT_11_1', 'MUT_15', 'WT_1']])
    bgrc = snps_readcount['WT_19'][s][alti]
    bgaf = snps_readcount['WT_19'][s][alti]/sum(snps_readcount['WT_19'][s])
    if bgrc <= 0.5*altrc and bgaf <= 0.5*altaf:
        goodw1.append((s, altrc, altaf, bgrc, bgaf))


################################################################################
############ SNP Identification when all eight samples are present #############
################################################################################
# Select SNPs in each sample of the eight sample:
readdepth = {'WT_1' : (70, 220),
             'wt7'  : (60, 160),
             'wt18' : (60, 160),
             'WT_19': (70, 250),
             'mut4' : (60, 150),
             'MUT_11_1' : (70, 250),
             'mut11_2': (60, 150),
             'MUT_15': (80, 230)}
samples = list(readdepth.keys())

## Test just the four new samples
### Checked positions that are present in either the mutant or the wt branches
### Result: No position in WT or MUT that have more than 20 reads supporting de novo mutation. So, no shared mutation between these branches.

## Testing all samples together
locwt1 = getlocs('WT_1/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=3)
locwt7 = getlocs('wt7/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=3)
locwt18 = getlocs('wt18/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=3)
locwt19 = getlocs('WT_19/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=3)
locmut4 = getlocs('mut4/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=3)
locmut11_1 = getlocs('MUT_11_1/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=3)
locmut11_2 = getlocs('mut11_2/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=3)
locmut15 = getlocs('MUT_15/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=3)

locwt1hi = getlocs('WT_1/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=20)
locwt7hi = getlocs('wt7/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=20)
locwt18hi = getlocs('wt18/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=20)
locwt19hi = getlocs('WT_19/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=20)
locmut4hi = getlocs('mut4/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=20)
locmut11_1hi = getlocs('MUT_11_1/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=20)
locmut11_2hi = getlocs('mut11_2/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=20)
locmut15hi = getlocs('MUT_15/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt', n=20)


locmthi = filterbg([locmut4hi, locmut11_1hi, locmut11_2hi, locmut15hi], [locwt1hi, locwt7hi, locwt18hi, locwt19hi], bgtype='any')
plot_selected_pos(locmthi, "tmp_igv.bat", "tmp/locmthi/", M=100, HEIGHT=200)
locmthi2 = filterbg([locmut4hi, locmut11_1hi, locmut11_2hi, locmut15hi], [locwt1, locwt7, locwt18, locwt19], bgtype='all')
#### Only high-confidence mutation is a two bp deletion at CUR4G_680602

#### Bootstrapping to see whether the mutation could be present by chance
loctst = filterbg([locwt1hi, locwt18hi, locmut4hi, locmut11_1hi], [locwt7hi, locwt19hi, locmut11_2hi, locmut15hi],  bgtype='any')
plot_selected_pos(loctst, "tmp_igv.bat", "tmp/loctst/", M=100, HEIGHT=200)
loctst = filterbg([locwt1hi, locwt18hi, locmut4hi, locmut11_1hi], [locwt7, locwt19, locmut11_2, locmut15],  bgtype='all')
loctst2 = filterbg([locwt7hi, locwt19hi, locmut11_2hi, locmut15hi], [locwt1hi, locwt18hi, locmut4hi, locmut11_1hi],  bgtype='any')
plot_selected_pos(loctst2, "tmp_igv.bat", "tmp/loctst2/", M=100, HEIGHT=200)
loctst2 = filterbg([locwt7hi, locwt19hi, locmut11_2hi, locmut15hi], [locwt1, locwt18, locmut4, locmut11_1], bgtype='all')
#### No mutation found in any config ==> the mutation in MUT sample is probably not random

#### Testing for mutation in WT samples absent in MUT sample
locwthi = filterbg([locwt1hi, locwt7hi, locwt18hi, locwt19hi], [locmut4hi, locmut11_1hi, locmut11_2hi, locmut15hi], bgtype='any')
plot_selected_pos(locwthi, "tmp_igv.bat", "tmp/locwthi/", M=75, HEIGHT=200)
locwthi2 = filterbg([locwt1hi, locwt7hi, locwt18hi, locwt19hi], [locmut4, locmut11_1, locmut11_2, locmut15], bgtype='all')
#### No mutation found. Probably this also mean that there is not any gene-coversion in the MUT samples.

#### Test mutations between branches that are closer to each other (that is diverged more recently)
##### Test-1: MUT11.1 and MUT11.2
l1_1 = filterbg([locmut11_1hi], [locwt1hi, locwt7hi, locwt18hi, locwt19hi, locmut4hi, locmut11_2hi, locmut15hi], bgtype='any')
l1_2 = filterbg([locmut11_1hi], [locwt1, locwt7, locwt18, locwt19, locmut4, locmut11_2, locmut15], bgtype='any')

l2_1 = filterbg([locmut11_2hi], [locwt1hi, locwt7hi, locwt18hi, locwt19hi, locmut4hi, locmut11_1hi, locmut15hi], bgtype='any')
l2_2 = filterbg([locmut11_2hi], [locwt1, locwt7, locwt18, locwt19, locmut4, locmut11_1, locmut15], bgtype='any')

l12_1 = filterbg([locmut11_1hi, locmut11_2hi], [locwt1hi, locwt7hi, locwt18hi, locwt19hi, locmut4hi, locmut15hi], bgtype='any')
l12_2 = filterbg([locmut11_1hi, locmut11_2hi], [locwt1, locwt7, locwt18, locwt19, locmut4, locmut15], bgtype='all')
l12_3 = filterbg([locmut11_1hi, locmut11_2hi], [locwt1, locwt7, locwt18, locwt19, locmut4, locmut15], bgtype='any')
plot_selected_pos(l12_2, "tmp_igv.bat", "tmp/l12_2/", M=75, HEIGHT=200)

##### Test-2: MUT4 and MUT15
l1_1 = filterbg([locmut4hi], [locwt1hi, locwt7hi, locwt18hi, locwt19hi, locmut11_1hi, locmut11_2hi, locmut15hi], bgtype='any')
l1_2 = filterbg([locmut4hi], [locwt1, locwt7, locwt18, locwt19, locmut11_1, locmut11_2, locmut15], bgtype='any')

l2_1 = filterbg([locmut15hi], [locwt1hi, locwt7hi, locwt18hi, locwt19hi, locmut4hi, locmut11_1hi, locmut11_2hi], bgtype='any')
l2_2 = filterbg([locmut15hi], [locwt1, locwt7, locwt18, locwt19, locmut4, locmut11_1, locmut11_2], bgtype='any')

l12_1 = filterbg([locmut4hi, locmut15hi], [locwt1hi, locwt7hi, locwt18hi, locwt19hi, locmut11_1hi, locmut11_2hi], bgtype='any')
l12_2 = filterbg([locmut4hi, locmut15hi], [locwt1, locwt7, locwt18, locwt19, locmut11_1, locmut11_2], bgtype='all')
l12_3 = filterbg([locmut4hi, locmut15hi], [locwt1, locwt7, locwt18, locwt19, locmut11_1, locmut11_2], bgtype='any')
plot_selected_pos(l12_2, "tmp_igv.bat", "tmp/l12_2/", M=75, HEIGHT=200)


##### Test-3: Testing pairs
p1 = filterbg([locmut4hi, locmut11_1hi], [locwt1, locwt7, locwt18, locwt19, locmut15, locmut11_2], bgtype='all')
p2 = filterbg([locmut15hi, locmut11_1hi], [locwt1, locwt7, locwt18, locwt19, locmut4, locmut11_2], bgtype='all')

p3 = filterbg([locmut4hi, locmut11_2hi], [locwt1, locwt7, locwt18, locwt19, locmut15, locmut11_1], bgtype='all')
p4 = filterbg([locmut15hi, locmut11_2hi], [locwt1, locwt7, locwt18, locwt19, locmut4, locmut11_1], bgtype='all')

p2filt = filterclosepos(p2)
p3filt = filterclosepos(p3)

p2filt =  filterbg([p2filt], [locwt1, locwt19], bgtype='all')
p3filt =  filterbg([p3filt], [locwt7, locwt18], bgtype='all')

p2filt =  filterbg([p2], [locwt1, locwt19], bgtype='all')
p3filt =  filterbg([p3], [locwt7, locwt18], bgtype='all')

################################################################################
######################### Indel Identification #################################
################################################################################

candidate_readcount = OrderedDict()
for sample in SAMPLES:
    print(sample)
    sample_dict = defaultdict()
    with open('{}/{}/multi_cell_{}_bt2_indel_candidate.sorted.filtered.bed'.format(INDIR, sample, sample), 'r') as fin:
        for line in fin:
            line = line.strip().split()
            sample_dict[line[0] + '_' + line[2]] = (line[3], line[4], int(line[5]), float(line[6]))
    candidate_readcount[sample] = sample_dict
candidate_readcount = OrderedDict(candidate_readcount)

loc = defaultdict(deque)
for sample in SAMPLES:
    for k, v in candidate_readcount[sample].items():
        loc[k].append(v[1])
for k in loc.keys():
    loc[k] = set(loc[k])


def get_readcounts(fin, loc):
    keys = set(loc.keys())
    sample_dict = {}
    with open(fin, 'r') as fin:
        for line in tqdm(fin):
            # print(line)
            # break
            if '+' not in line:
                if '-' not in line:
                    continue
            line = line.strip().split()
            if line[0] + '_' + line[1] in keys:
                sample_dict[line[0] + '_' + line[1]] = {line[i]: int(line[i+1]) for i in range(9, len(line), 2) if line[i] in loc[line[0] + '_' + line[1]]}
                sample_dict[line[0] + '_' + line[1]]['RC'] = int(line[3])
    return sample_dict

fins = ['{}{}/bam_read_counts_b30_q10.bt2.txt'.format(INDIR, sample) for sample in SAMPLES]
with Pool(processes=4) as pool:
    snps_readcount = pool.map(partial(get_readcounts, loc=loc), fins)
snps_readcount = dict(zip(SAMPLES, snps_readcount))

# Get indels present in both MUT or both WT
mtloc = deque()
for l, v in candidate_readcount['MUT_11_1'].items():
    if v[3] < 0.1: continue
    if v[2] < 5: continue
    try:
        if snps_readcount['MUT_15'][l][v[1]] < 5: continue
        if snps_readcount['MUT_15'][l][v[1]]/snps_readcount['MUT_15'][l]['RC'] < 0.1: continue
    except KeyError:
        continue
    mtloc.append(l)

badloc = deque()
for l in mtloc:
    v = candidate_readcount['MUT_11_1'][l][1]
    altrc = min(snps_readcount['MUT_11_1'][l][v], snps_readcount['MUT_15'][l][v])
    try:
        w1 = int(snps_readcount['WT_1'][l][v])
    except KeyError:
        w1 = 0
    try:
        w2 = int(snps_readcount['WT_19'][l][v])
    except KeyError:
        w2 = 0
    if w1 > 0.2*altrc and w2 > 0.2*altrc:       # using AND operator instead of OR. Should change back later
        badloc.append(l)
mtloc2 = set(mtloc) - set(badloc)
mtloc2 = [i.rsplit("_", 1) for i in mtloc2]

# Filter mutations that have another indel within 100bp
chrs = set([i[0] for i in mtloc2])
chrdic = defaultdict(deque)
for v in mtloc2:
    chrdic[v[0]].append(int(v[1]))
for k, v in chrdic.items():
    vs = sorted(v)
    bp = deque()
    for i in range(1, len(vs)):
        if vs[i] - vs[i-1] < 100:
            bp.append(vs[i-1])
            bp.append(vs[i])
    bp = set(bp)
    gp = sorted(list(set(vs) - set(bp)))
    chrdic[k] = gp
mtloc2 = [k+"_"+str(p) for k, v in chrdic.items() for p in v]

# Write de-novo indels
with open(CWD + "denovo_indels.txt", 'w') as fout:
    fout.write("chr{t}pos{t}ref{t}alt{t}wt1_rc{t}wt1_af{t}wt19_rc{t}wt19_af{t}mut11_rc{t}mut11_af{t}mut15_rc{t}mut15_af".format(t="\t") + "\n")
    for l in mtloc2:
        refi = BASE_DICT[candidate_readcount['MUT_11_1'][l][0]]
        alt = candidate_readcount['MUT_11_1'][l][1]
        outstr = []
        outstr += l.rsplit("_", 1)
        outstr += list(candidate_readcount['MUT_11_1'][l][:2])
        try:
            outstr += [snps_readcount['WT_1'][l][alt], str(round(int(snps_readcount['WT_1'][l][alt])/int(snps_readcount['WT_1'][l]['RC']), 4))]
        except KeyError:
            outstr += ['0', '0']
        try:
            outstr += [snps_readcount['WT_19'][l][alt], str(round(int(snps_readcount['WT_19'][l][alt])/int(snps_readcount['WT_19'][l]['RC']), 4))]
        except KeyError:
            outstr += ['0', '0']
        try:
            outstr += [snps_readcount['MUT_11_1'][l][alt], str(round(int(snps_readcount['MUT_11_1'][l][alt])/int(snps_readcount['MUT_11_1'][l]['RC']), 4))]
        except KeyError:
            outstr += ['0', '0']
        try:
            outstr += [snps_readcount['MUT_15'][l][alt], str(round(int(snps_readcount['MUT_15'][l][alt])/int(snps_readcount['MUT_15'][l]['RC']), 4))]
        except KeyError:
            outstr += ['0', '0']
        fout.write("\t".join(outstr) + "\n")

wtloc = deque()
for l, v in candidate_readcount['WT_1'].items():
    if v[2] < 5: continue
    try:
        if int(snps_readcount['WT_19'][l][v[1]]) < 5: continue
    except KeyError:
        continue
    wtloc.append(l)

badloc = deque()
for l in wtloc:
    v = candidate_readcount['WT_1'][l][1]
    altrc = min(int(snps_readcount['WT_1'][l][v]), int(snps_readcount['WT_19'][l][v]))
    try:
        w1 = int(snps_readcount['MUT_11_1'][l][v])
    except KeyError:
        w1 = 0
    try:
        w2 = int(snps_readcount['MUT_15'][l][v])
    except KeyError:
        w2 = 0
    if w1 > 0.2*altrc or w2 > 0.2*altrc:
        badloc.append(l)
wtloc2 = set(wtloc) - set(badloc)
wtloc2 = [i.rsplit("_", 1) for i in wtloc2]

# Filter mutations that have another indel within 100bp
chrs = set([i[0] for i in wtloc2])
chrdic = defaultdict(deque)
for v in wtloc2:
    chrdic[v[0]].append(int(v[1]))
for k, v in chrdic.items():
    vs = sorted(v)
    bp = deque()
    for i in range(1, len(vs)):
        if vs[i] - vs[i-1] < 100:
            bp.append(vs[i-1])
            bp.append(vs[i])
    bp = set(bp)
    gp = sorted(list(set(vs) - set(bp)))
    chrdic[k] = gp
wtloc2 = [k+"_"+str(p) for k, v in chrdic.items() for p in v]

# Write indels present in the WT but lost in the MUTs
with open(CWD + "indel_lose.txt", 'w') as fout:
    fout.write("chr{t}pos{t}ref{t}alt{t}wt1_rc{t}wt1_af{t}wt19_rc{t}wt19_af{t}mut11_rc{t}mut11_af{t}mut15_rc{t}mut15_af".format(t="\t") + "\n")
    for l in wtloc2:
        refi = BASE_DICT[candidate_readcount['WT_1'][l][0]]
        alt = candidate_readcount['WT_1'][l][1]
        outstr = []
        outstr += l.rsplit("_", 1)
        outstr += list(candidate_readcount['WT_1'][l][:2])
        try:
            outstr += [snps_readcount['WT_1'][l][alt], str(round(int(snps_readcount['WT_1'][l][alt])/int(snps_readcount['WT_1'][l]['RC']), 4))]
        except KeyError:
            outstr += ['0', '0']
        try:
            outstr += [snps_readcount['WT_19'][l][alt], str(round(int(snps_readcount['WT_19'][l][alt])/int(snps_readcount['WT_19'][l]['RC']), 4))]
        except KeyError:
            outstr += ['0', '0']
        try:
            outstr += [snps_readcount['MUT_11_1'][l][alt], str(round(int(snps_readcount['MUT_11_1'][l][alt])/int(snps_readcount['MUT_11_1'][l]['RC']), 4))]
        except KeyError:
            outstr += ['0', '0']
        try:
            outstr += [snps_readcount['MUT_15'][l][alt], str(round(int(snps_readcount['MUT_15'][l][alt])/int(snps_readcount['MUT_15'][l]['RC']), 4))]
        except KeyError:
            outstr += ['0', '0']
        fout.write("\t".join(outstr) + "\n")

################################################################################
########################## TE Indels Identification ############################
################################################################################
# Analysing TE identified using TEPID
CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/'
SAMPLES = ("MUT_11_1", "MUT_15", "WT_1", "WT_19")

## Deletions
deldf = pd.DataFrame()
for s in SAMPLES:
    df = pd.read_table(CWD+"{s}/deletions_{s}.bed".format(s=s), header=None)
    df['s'] = s
    deldf = pd.concat([deldf, df])
upsetdf = from_contents({'WT_1': pd.unique(deldf.loc[deldf['s'] == 'WT_1', 4]),
                         'WT_19': pd.unique(deldf.loc[deldf['s'] == 'WT_19', 4]),
                         'MUT_11_1': pd.unique(deldf.loc[deldf['s'] == 'MUT_11_1', 4]),
                         'MUT_15': pd.unique(deldf.loc[deldf['s'] == 'MUT_15', 4])
                         })
UpSet(upsetdf, subset_size='count').plot()
fig = plt.figure(figsize=[8, 8])
UpSet(upsetdf, subset_size='count').plot(fig=fig)
plt.show()
plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/tepid_TE_upset_deletions.pdf')
plt.close()
wtdel = deque()
mutdel = deque()
for g in deldf.groupby(4):
    if len(g[1]) != 2: continue
        s = list(g[1]['s'])
        if 'WT_1' in s and 'WT_19' in s: wtdel.append(g[0])
        if 'MUT_11_1' in s and 'MUT_15' in s: mutdel.append(g[0])
mutdeldf = deldf.loc[deldf[4].isin(mutdel)].sort_values([4])
wtdeldf = deldf.loc[deldf[4].isin(wtdel)].sort_values([4])
mutdeldf.to_csv("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/mut_del_tepid_TE.txt", header=False, index=False, sep='\t')
wtdeldf.to_csv("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/wt_del_tepid_TE.txt", header=False, index=False, sep='\t')
tepid_te_del_igv_batch(CWD+'mut_del_tepid_TE.txt', CWD+'mut_del_tepid_TE.igv.batch', CWD+'mut_del_tepid_TE_snapshots')

## Insertions
insdf = pd.DataFrame()
wt1 = bdt.BedTool(CWD+"{s}/insertions_{s}.bed".format(s='WT_1'))
wt19 = bdt.BedTool(CWD+"{s}/insertions_{s}.bed".format(s='WT_19'))
mut11 = bdt.BedTool(CWD+"{s}/insertions_{s}.bed".format(s='MUT_11_1'))
mut15 = bdt.BedTool(CWD+"{s}/insertions_{s}.bed".format(s='MUT_15'))
BPDIFF = 15
names = ['chr', 's', 'e', 'cin', 'sin', 'ein', 'name', 'id', 'chr2', 's2', 'e2', 'cin2', 'sin2', 'ein2', 'name2', 'id2', 'l']
mut  = mut11.intersect(mut15, wao=True).to_dataframe(names=names)
m1w1 = mut11.intersect(wt1, wao=True).to_dataframe(names=names)
m1w2 = mut11.intersect(wt19, wao=True).to_dataframe(names=names)
m2w1 = mut15.intersect(wt1, wao=True).to_dataframe(names=names)
m2w2 = mut15.intersect(wt19, wao=True).to_dataframe(names=names)

w1df = wt1.to_dataframe(names=['chr', 's', 'e', 'cin', 'sin', 'ein', 'name', 'id'])
w19df = wt19.to_dataframe(names=['chr', 's', 'e', 'cin', 'sin', 'ein', 'name', 'id'])


def chkovrlap(m , df, t):
    if len(df.loc[(df['chr']==m[0]) &
           ((df['s'] - m[2]) < t) &
           ((m[1] - df['e']) < t)]) > 0:
        return False
    if len(df.loc[(df['chr']==m[8]) &
           ((df['s'] - m[10]) < t) &
           ((m[9] - df['e']) < t)]) > 0:
        return False
    return True

goodmt = deque()
test = deque()
for m in mut.itertuples(index=False):
    if not any([p in m[14].split(',') for p in m[6].split(',')]): continue
    if m[1]-m[9] > BPDIFF and m[2]-m[10] > BPDIFF: continue
    if m[9]-m[1] > BPDIFF and m[10]-m[2] > BPDIFF : continue
    if not chkovrlap(m, w1df, 1000): continue
    if not chkovrlap(m, w19df, 1000): continue
    goodmt.append(m)

mutinsdf = pd.DataFrame(goodmt)
mutinsdf.to_csv("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/mut_ins_tepid_TE.txt", header=False, index=False, sep='\t')
tepid_te_ins_igv_batch(CWD+'mut_ins_tepid_TE.txt', CWD+'mut_ins_tepid_TE.igv.batch', CWD+'mut_ins_tepid_TE_snapshots')

### Checking to see if there are branch specific TE insertions
BPDIFF = 15
names = ['chr', 's', 'e', 'cin', 'sin', 'ein', 'name', 'id', 'chr2', 's2', 'e2', 'cin2', 'sin2', 'ein2', 'name2', 'id2', 'l']
mut11m = bdt.BedTool([Interval(f.chrom, max(0, f.start-500), f.end+500) for f in mut11]).merge()
mut15m = bdt.BedTool([Interval(f.chrom, max(0, f.start-500), f.end+500) for f in mut15])
wt1m = bdt.BedTool([Interval(f.chrom, max(0, f.start-500), f.end+500) for f in wt1])
wt19m = bdt.BedTool([Interval(f.chrom, max(0, f.start-500), f.end+500) for f in wt19])

bg = mut15m.cat(wt1m, postmerge=True).cat(wt19m, postmerge=True)
mt11sp = mut11m.intersect(bg, v=True)
mt11spori = mut11.intersect(mt11sp, u=True)
mt11sp.saveas(CWD+'/MUT_11_1/MUT_11_1_specific_TE_ins.bed')
tepid_te_ins_bs_igv_batch(CWD+'/MUT_11_1/MUT_11_1_specific_TE_ins.bed', CWD+'/MUT_11_1/igv_ins.bat', CWD+'/MUT_11_1/ins_snap')
### On manual checking, NO MUT_11_1 specific TE insertion was identified. So, other samples were not checked.

# Analysing TE identified using TEfinder
## Get TEs that PASS the filtering
def filtered_grange(bed, m=100):
    """
    also adds margin (increase size) on both ends of the range
    """
    fin=pr.read_bed(bed)
    fin.columns = ['Chromosome', 'Start', 'End', 'Name', 'read_cnt', 'bias', 'info']
    # finf = fin[pd.Series('PASS' in i for i in fin.info)]
    finf = fin
    finf.Start = finf.Start - m
    finf.End = finf.End + m
    return finf

CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/tefinder/'
mut11 = filtered_grange(CWD+'MUT_11_1/TEinsertions.sorted.bed', m=500)
mut15 = filtered_grange(CWD+'MUT_15/TEinsertions.sorted.bed', m=500)
wt1 = filtered_grange(CWD+'WT_1/TEinsertions.sorted.bed', m=500)
wt19 = filtered_grange(CWD+'WT_19/TEinsertions.sorted.bed', m=500)

mt11f = mut11[pd.Series('PASS' in i for i in mut11.info)]
mt15f = mut15[pd.Series('PASS' in i for i in mut15.info)]
wt1f = wt1[pd.Series('PASS' in i for i in wt1.info)]
wt19f = wt19[pd.Series('PASS' in i for i in wt19.info)]

mt11o = mt11f.overlap(mt15f)
mt15o = mt15f.overlap(mt11f)
mt = pr.concat([mt11o, mt15o]).merge()

wt1o = wt1f.overlap(wt19f)
wt19o = wt19f.overlap(wt1f)
wt = pr.concat([wt1o, wt19o]).merge()

mtsp = mt.overlap(pr.concat([wt1, wt19]).merge(), invert=True)
wtsp = wt.overlap(pr.concat([mut11, mut15]).merge(), invert=True)
### this results in 2 and 8 mut/wt specific TE insertions. But, in IGV all branches have TE insertions reads

### Checking to see if there are branch specific TE insertions
mt11sp = mut11.overlap(pr.concat([mut15, wt1, wt19]).merge(), invert=True)
mt11spf = mt11sp[pd.Series((mt11sp.read_cnt > 50) & (mt11sp.read_cnt < 300))].sort(by='read_cnt')
pos = mt11spf.Chromosome.astype(str) + ":" + mt11spf.Start.astype(str) + "-" + mt11spf.End.astype(str)
### On manual checking, NO MUT_11_1 specific TE insertion was identified. So, other samples were not checked.

################################################################################
########################## Pathogen infection test #############################
################################################################################

# Generate species vector for each branch
fname = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/patho_inf/{}_megahit_assembly/ntall.oblast'
spcnt = {}
for s in ('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15'):
    cnts = defaultdict(int)
    df = pd.read_table(fname.format(s), header=None)
    for grp in df.groupby([0]):
        sps = set(["_".join(i.split()[:2]) for i in grp[1][16]])
        for sp in sps:
            cnts[sp] += 1
    spcnt[s] = cnts

allsp = set(list(spcnt['WT_1'])+list(spcnt['WT_19'])+list(spcnt['MUT_11_1'])+list(spcnt['MUT_15'])
allspcnt = {s:[spcnt[p][s] for p in ('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15')] for s in allsp}
varsp = {s:allspcnt[s] for s in allsp if allspcnt[s][0]+allspcnt[s][1] < 0.5*(allspcnt[s][2]+allspcnt[s][3]) or allspcnt[s][0]+allspcnt[s][1] > 2*(allspcnt[s][2]+allspcnt[s][3])}
varsphigh = {k:v for k,v in varsp.items() if any([True if i > 10 else False for i in v])}
### print(varsphigh)
### {'Primula_capitata': [8, 11, 3, 6],
###  'Ixonanthes_chinensis': [8, 17, 6, 5],
###  'Pontederia_cordata': [7, 13, 5, 4],
###  'Maranta_leuconeura': [8, 13, 4, 6],
###  'Homo_sapiens': [9, 14, 5, 5]}
