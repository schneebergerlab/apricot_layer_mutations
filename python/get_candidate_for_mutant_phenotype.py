import os
from collections import defaultdict, OrderedDict, deque
from tqdm import tqdm
from multiprocessing import Pool
from functools import partial
import pandas as pd
from subprocess import run
import sys
sys.path.insert(0, '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python')
from get_candidate_for_mutant_phenotype_funcs import getlinefrombamrc
from multiprocessing import set_start_method, Pool
# set_start_method('spawn')
from matplotlib import pyplot as plt
import numpy as np
from datetime import datetime

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

# positions which have either double or half allele-frequency in MUT vs WT
goodmt = deque()
for s in loc:
    try:
        wrc = snps_readcount['WT_1'][s][0] + snps_readcount['WT_19'][s][0]
        waf = wrc/(snps_readcount['WT_1'][s][1] + snps_readcount['WT_19'][s][1])
        mrc = snps_readcount['MUT_11_1'][s][0] + snps_readcount['MUT_15'][s][0]
        maf = mrc/(snps_readcount['MUT_11_1'][s][1] + snps_readcount['MUT_15'][s][1])
        if wrc < 0.5*mrc and waf < 0.5*maf:
            goodmt.append([s, wrc, waf, mrc, maf])
        elif wrc > 2*mrc and waf > 2*maf:
            goodmt.append([s, wrc, waf, mrc, maf])
    except ZeroDivisionError:
        pass

x0, x1 = 0, 1
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
for s in goodmt:
    if s[1] < 20 and s[3] < 20: continue
    # fig.add_artist(lines.Line2D([0, 1], (s[1], s[3])))
    ax.plot([0, 1], (s[1], s[3]))


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
def foo():
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
    cnt = 0
    while True:
        print(cnt, str(datetime.now()))
        cnt+=1
        if cnt == 5: break
        changed = False
        for k, v in files.items():
            for i in range(READSIZE):
                w1l = v.readline().strip().split()
                if w1l == '': break
                else:
                    x = list(map(int, w1l[4:8]))
                    posdata[(w1l[0], w1l[1], w1l[2])][k] = x
                    changed = True
        rem = deque()
        for k, v in posdata.items():
            if len(v) < 4:
                #TODO: FIX the below line
                if k[0] != w1l[0]: rem.append(k)
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
with open("TMP_goodpos.pickle", 'rb') as fin:
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

fig = plt.figure()
ax = fig.add_subplot(2,1,1)
ax.scatter([v[0] for v in posdata.values()], [v[2] for v in posdata.values()])
ax.set_xlim(-10, 200)
ax.set_ylim(-10, 200)
ax.set_title("Alternate allele read-count")
ax.set_xlabel("Merged WT")
ax.set_ylabel("Merged MUT")

ax = fig.add_subplot(2,1,2)
ax.scatter([v[0]/v[1] for v in posdata.values()], [v[2]/v[3] for v in posdata.values()])
ax.set_xlim(-0.05, 1)
ax.set_ylim(-0.05, 1)
ax.set_title("Alternate allele frequency")
ax.set_xlabel("Merged WT")
ax.set_ylabel("Merged MUT")

fig = plt.figure()
ax = fig.add_subplot(2,1,1)
ax.scatter([v[0] for v in garb.values()], [v[2] for v in garb.values()])
ax.set_xlim(-10, 50)
ax.set_ylim(-10, 200)
ax.set_title("Alternate allele read-count")
ax.set_xlabel("Merged WT")
ax.set_ylabel("Merged MUT")

ax = fig.add_subplot(2,1,2)
ax.scatter([v[0]/v[1] for v in garb.values()], [v[2]/v[3] for v in garb.values()])
ax.set_xlim(-0.05, 1)
ax.set_ylim(-0.05, 1)
ax.set_title("Alternate allele frequency")
ax.set_xlabel("Merged WT")
ax.set_ylabel("Merged MUT")



chrlen=OrderedDict()
with open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.chrlen.txt", "r") as fin:
    for l in fin:
        l = l.strip().split()
        chrlen[l[0]] = int(l[1])



cand = deque()
w1l, w2l, m1l, m2l = None, None, None, None
for c in chrlen.keys():
    rcarr = np.array(['']*(chrlen[c]+1), dtype=object)
    if w1l is not None: rcarr[int(w1l[1])] += ",".join(w1l[4:8])+';'
    for line in tqdm(w1):
        line = line.strip().split()
        if line[0] == c:
            rcarr[int(line[1])] += ",".join(line[4:8])+';'
        else:
            w1l = line
            break
    if w2l is not None: rcarr[int(w2l[1])] += ",".join(w2l[4:8])+';'
    for line in tqdm(w2):
        line = line.strip().split()
        if line[0] == c:
            rcarr[int(line[1])] += ",".join(line[4:8])+';'
        else:
            w2l = line
            break
    if m1l is not None: rcarr[int(m1l[1])] += ",".join(m1l[4:8])+';'
    for line in tqdm(m1):
        line = line.strip().split()
        if line[0] == c:
            rcarr[int(line[1])] += ",".join(line[4:8])+';'
        else:
            m1l = line
            break
    if m2l is not None: rcarr[int(m2l[1])] += ",".join(m2l[4:8])+';'
    for line in tqdm(m2):
        line = line.strip().split()
        if line[0] == c:
            rcarr[int(line[1])] += ",".join(line[4:8])+';'
        else:
            m2l = line
            break
    for p in rcarr:
        rcw1
    break


p1 = w1.readline().strip().split()
p2 = w2.readline().strip().split()
p3 = m1.readline().strip().split()
p4 = m2.readline().strip().split()
while True:
    if p1 == '' and p2 == '' and p3 == '' and p4 == '': break
    if p1[0] == p2[0] == p3[0] == p4[0]: c = p1[0]



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

