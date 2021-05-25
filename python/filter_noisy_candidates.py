import os
import pandas as pd
from subprocess import run
from collections import deque
from multiprocessing import Pool
from functools import partial

SAMPLES     = ("WT_1", "WT_19", "MUT_11_1", "MUT_15")
CWD         = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
T           = 100 # Threashold for neighbour

os.chdir(CWD)
for sample in SAMPLES:
    os.chdir(CWD+sample)
    with open(sample + "_good_candidates.txt", 'r') as fin:
        with open(sample + "_good_candidates_nonear.txt", 'w') as fout:
            cur = []
            curbad = False
            for line in fin:
                line = line.strip().split()
                if cur == []:
                    cur = line
                else:
                    if line[0] != cur[0]:
                        if not curbad:
                            fout.write('\t'.join(cur) + "\n")
                        cur = line
                        curbad = False
                    else:
                        if int(line[2]) - int(cur[2]) >= 100:
                            if not curbad:
                                fout.write('\t'.join(cur) + "\n")
                            cur = line
                            curbad = False
                        else:
                            cur = line
                            curbad = True
            if not curbad:
                fout.write('\t'.join(cur) + "\n")

os.chdir(CWD)
dfs = []
for s in SAMPLES:
    d = pd.read_table(CWD + s + "/" + s + "_good_candidates_nonear.txt", header=None, index_col=False)
    d[7] = s
    dfs.append(d)
df = pd.concat(dfs, sort=True)
df.sort_values([0, 1, 2], inplace=True)
df = df.loc[df[5] < 10]
df[1] = df[1] + 1
df.to_csv(CWD+'all_good_candidate.txt', index=False, header=False, sep=' ')

for sample in SAMPLES:
    os.chdir(CWD + sample)
    f = open("TMP.txt", 'w')
    run("bam-readcount -w 0 -b 0 -q 0 -l ../all_good_candidate.txt -f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta {}.sorted.bt2.bam ".format(sample).split(), stdout=f)
    f.close()
    with open("TMP.txt", 'r') as fin:
        with open('all_good_readcount.b0.q0.txt', 'w') as fout:
            for line in fin:
                line = line.strip().split()
                fout.write("\t".join(line[:4] + [i.split(':')[1] for i in line[5:11]])+"\n")

    f = open("TMP.txt", 'w')
    run("bam-readcount -w 0 -b 30 -q 40 -l ../all_good_candidate.txt -f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta {}.sorted.bt2.bam ".format(sample).split(), stdout=f)
    f.close()
    with open("TMP.txt", 'r') as fin:
        with open('all_good_readcount.b30.q40.txt', 'w') as fout:
            for line in fin:
                line = line.strip().split()
                fout.write("\t".join(line[:4] + [i.split(':')[1] for i in line[5:11]])+"\n")
    os.remove('TMP.txt')

# Filter positions for which other samples get >5 reads when no mapping quality filter is used
# OR which have <3 reads when mapping quality > 40
dfnomap = {}
dfhmap = {}
BASE_DICT= {'A': 4, 'C': 5, 'G': 6, 'T': 7}
for sample in SAMPLES:
    d = pd.read_table("{}/{}/all_good_readcount.b0.q0.txt".format(CWD, sample), header=None, engine='python', index_col=False)
    d[9] = sample
    dfnomap[sample] = d

    d = pd.read_table("{}/{}/all_good_readcount.b30.q40.txt".format(CWD, sample), header=None, engine='python', index_col=False)
    d[9] = sample
    dfhmap[sample] = d

isbad = deque()
cnt1 = 0
cnt2 = 0
for i in range(df.shape[0]):
    s = df.iat[i, 7]
    # if dfhmap[s].iat[i, BASE_DICT[df.iat[i, 4]]] < 5:
    if dfhmap[s].iat[i, BASE_DICT[df.iat[i, 4]]] < 3:
        isbad.append(True)
        cnt1 += 1
        continue
    ac = 0
    for sample in SAMPLES:
        if sample == s: continue
        ac += dfnomap[sample].iat[i, BASE_DICT[df.iat[i, 4]]]
    if ac >= 5:
        isbad.append(True)
        cnt2 += 1
        continue
    isbad.append(False)

baddf = df.iloc[[i for i in range(len(isbad)) if isbad[i] is True]].copy()
gooddf = df.iloc[[i for i in range(len(isbad)) if isbad[i] is False]].copy()
gooddf.to_csv(CWD+'all_good_candidate.qualfilter.txt', index=False, header=False, sep=' ')
gooddf[1] -= 1
gooddf.to_csv(CWD+'all_good_candidate.qualfilter.bed', index=False, header=False, sep=' ')

def samplempileup(sample, CWD):
    print(sample)
    os.chdir(CWD + sample)
    # run("/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools view -b -o all_good_candidate.qualfilter.bam -L ../all_good_candidate.qualfilter.bed {}.sorted.bt2.bam".format(sample).split())
    # with open('all_good_candidate.qualfilter.baq.bam', 'wb') as f:
    #     run("/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools calmd -Abr all_good_candidate.qualfilter.bam /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta".format(sample).split(), stdout=f)
    # run("/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools index all_good_candidate.qualfilter.baq.bam".split())
    os.remove('all_good_candidate.qualfilter.baq.mpileup')
    with open('all_good_candidate.qualfilter.baq.mpileup', 'a') as fout:
        for row in gooddf.itertuples(index=False):
            run("/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools mpileup -aa -A -f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta -q 40 -Q 13 -r {}:{}-{} all_good_candidate.qualfilter.baq.bam".format(row[0], row[1], row[2]).split(), stdout=fout)


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


isgood = [False]*gooddf.shape[0]
for sample in SAMPLES:
    os.chdir(CWD + sample)
    mpileup = pd.read_table('all_good_candidate.qualfilter.baq.mpileup', header=None, index_col=False)
    for i in range(gooddf.shape[0]):
        if gooddf.iat[i, 7] == sample:
            _, bases = getbases(mpileup.iat[i, 4])
            if len([c for c in bases if c in [gooddf.iat[i, 4].upper(), gooddf.iat[i, 4].lower()]]) >= 5:
                isgood[i] = True
                continue

low_cov_mutants = gooddf.iloc[isgood].copy()
low_cov_mutants.to_csv(CWD+"/low_cov_mutants.txt", header=False, index=False)
# The selected 69 candidates were manually curated to get high confidence low coverage somatic mutations