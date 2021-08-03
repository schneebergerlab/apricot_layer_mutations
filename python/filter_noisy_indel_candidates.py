import os
import pandas as pd
from subprocess import run
from collections import deque, Counter, defaultdict
from multiprocessing import Pool
from functools import partial
import pysam
import sys
sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/software/hometools/')
from myUsefulFunctions import cggenlen
from tqdm import tqdm
import pickle
# sys.path.insert(0, '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/')
# from filter_noisy_candidates import getbases


SAMPLES     = ("WT_1", "WT_19", "MUT_11_1", "MUT_15")
CWD         = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
T           = 100 # Threashold for neighbour

os.chdir(CWD)
for sample in SAMPLES:
    os.chdir(CWD+sample)
    with open(sample + "_indel_good_candidates.txt", 'r') as fin:
        with open(sample + "_indel_good_candidates_nonear.txt", 'w') as fout:
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
                        if int(line[2]) - int(cur[2]) >= T:
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
    d = pd.read_table(CWD + s + "/" + s + "_indel_good_candidates_nonear.txt", header=None, index_col=False)
    d[7] = s
    dfs.append(d)
df = pd.concat(dfs, sort=True)
df.sort_values([0, 1, 2], inplace=True)
df = df.loc[df[5] < 10]
df[1] = df[1] + 1
df.to_csv(CWD+'all_good_candidate_indels.txt', index=False, header=False, sep=' ')

garb = df.copy()
garb.reset_index(inplace=True, drop=True)
delin = ['-' in g for g in garb[4]]
insin = ['+' in g for g in garb[4]]
garb.loc[delin, 1] = garb.loc[delin, 1].astype(int) - 1
garb.loc[insin, 2] = garb.loc[insin, 2].astype(int) + 1
garb.to_csv(CWD+'indel_loci_for_readcount.txt', index=False, header=False, sep=' ')

# awk '{printf $1" "$2" "$3" "$4; for(i=6;i<=10;i++) {n1=split($i,a,":"); printf " "a[2]};  for(i=11;i<=NF;i++) {n1=split($i,a,":"); printf " "a[1]" "a[2]}; printf "\n"}' > $out



def readcount(sample, CWD):
    os.chdir(CWD + sample)
    f = open("TMP.txt", 'w')
    run("bam-readcount -w 0 -b 0 -q 0 -l ../indel_loci_for_readcount.txt -f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta {}.sorted.bt2.bam ".format(sample).split(), stdout=f)
    f.close()
    with open("TMP.txt", 'r') as fin:
        with open('all_good_readcount_indels.b0.q0.txt', 'w') as fout:
            for line in fin:
                line = line.strip().split()
                fout.write("\t".join(line[:4] + [i.split(':')[1] for i in line[5:10]] + [j for i in line[10:] for j in i.split(':')[:2] ]) + "\n")

    f = open("TMP.txt", 'w')
    run("bam-readcount -w 0 -b 30 -q 40 -l ../indel_loci_for_readcount.txt -f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta {}.sorted.bt2.bam ".format(sample).split(), stdout=f)
    f.close()
    with open("TMP.txt", 'r') as fin:
        with open('all_good_readcount_indels.b30.q40.txt', 'w') as fout:
            for line in fin:
                line = line.strip().split()
                fout.write("\t".join(line[:4] + [i.split(':')[1] for i in line[5:10]] + [j for i in line[10:] for j in i.split(':')[:2]]) + "\n")
    os.remove('TMP.txt')


with Pool(processes=2) as pool:
    pool.map(partial(readcount, CWD=CWD), SAMPLES)


# Filter positions for which other samples get >5 reads when no mapping quality filter is used
# OR which have <3 reads when mapping quality > 40
def readindelcounts(fin):
    sindel = dict()
    with open(fin, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            if len(line) < 10: continue
            sindel[(line[0], int(line[1]))] = {line[i]: int(line[i+1]) for i in range(9, len(line), 2)}
    return sindel

allindels = dict(zip(zip(df[0], df[2]), df[4]))
lsindel = {sample: readindelcounts(CWD + sample + "/all_good_readcount_indels.b0.q0.txt") for sample in SAMPLES}
hsindel = {sample: readindelcounts(CWD + sample + "/all_good_readcount_indels.b30.q40.txt") for sample in SAMPLES}

goodind = deque()
for row in df.itertuples(index=False):
    try: hc = hsindel[row[7]][(row[0], row[1])][row[4]]
    except KeyError: hc = 0
    if hc < 3: continue
    tmp = row[2]-1 if '-' in row[4] else row[2]+1
    loci = set([(sample, pos) for sample in SAMPLES for pos in [row[2], tmp]])
    lr = 0
    for l in loci:
        try: lr += sum(lsindel[l[0]][(row[0], l[1])].values())
        except KeyError: pass
    if lr - lsindel[row[7]][(row[0], row[2])][row[4]] >= 5: continue
    goodind.append(row)

gooddf = pd.DataFrame(goodind, index=None, columns=range(8))
gooddf[1] = gooddf[1] - 1
gooddf[list(range(8))].to_csv(CWD+'all_good_candidate_indels.qualfilter.bed', index=False, header=False, sep=' ')
gooddf[1] = gooddf[1] + 1


def samplempileup(sample, CWD):
    print(sample)
    os.chdir(CWD + sample)
    # run("/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools view -b -o all_good_candidate_indels.qualfilter.bam -L ../all_good_candidate_indels.qualfilter.bed {}.sorted.bt2.bam".format(sample).split())
    # with open('all_good_candidate_indels.qualfilter.baq.bam', 'wb') as f:
    #     run("/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools calmd -Abr all_good_candidate_indels.qualfilter.bam /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta".format(sample).split(), stdout=f)
    # run("/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools index all_good_candidate_indels.qualfilter.baq.bam".split())

    try:
        os.remove('all_good_candidate_indels.qualfilter.baq.mpileup')
    except FileNotFoundError:
        pass
    with open('all_good_candidate_indels.qualfilter.baq.mpileup', 'w') as fout:
        # for row in df.itertuples(index=False):
        run("/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools mpileup -f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta -l ../TMP.bed -aa -A -q 40 -Q 0 --rf 3 --reverse-del --output-QNAME -O --output-extra RG,FLAG all_good_candidate_indels.qualfilter.baq.bam".split(), stdout=fout)
#end

with Pool(processes=2) as pool:
    # Get TMP.bed to allow selection of position before deletion as well
    with open(CWD+'TMP.bed', 'w') as fout:
        with open(CWD + 'all_good_candidate_indels.qualfilter.bed', 'r') as fin:
            for line in fin:
                line = line.strip().split()
                line[1] = str(int(line[1]) - 1)
                fout.write(' '.join(line) + '\n')
    out = pool.map(partial(samplempileup, CWD=CWD), SAMPLES)
    os.remove(CWD+'TMP.bed')

def goodmapread(bamfin):
    bam = pysam.AlignmentFile(bamfin, 'rb')
    sreadid = deque()
    references = bam.references
    for r in tqdm(bam):
        if r.is_unmapped: continue
        if r.mate_is_unmapped: continue
        if r.is_reverse == r.mate_is_reverse: continue
        if r.reference_name != references[r.mrnm]: continue
        if abs(r.reference_start - r.mpos) > 1500: continue
        if r.is_reverse:
            if r.mpos > r.reference_start: continue
        else:
            if r.mpos < r.reference_start: continue
        try: mrlen = cggenlen(r.opt('MC'), 'r')
        except KeyError: continue
        if r.is_reverse:
            if r.mpos + mrlen > r.reference_end: continue
        else:
            if r.mpos + mrlen < r.reference_end: continue
        sreadid.append([r.qname, r.reference_name, r.reference_start, r.reference_end, r.is_reverse, r.is_read1, r.mpos, r.mpos + mrlen])
    return sreadid
#end

def getindels(l):
    indeldata = deque()
    pos = 0
    ind = -1
    for i in range(len(l)):
        if i < pos: continue
        if l[i] == '$': continue
        if l[i] == '^':
            ind -= 1
            continue
        ind += 1
        # print(l[i], ind)
        if l[i] in ['-', '+']:
            ind -= 1
            for j in range(i+1, len(l)):
                if l[j] not in nums: break
            ilen = int(l[(i+1):j])
            iseq = l[j: (j+ilen)]
            index = ind
            pos = j + ilen
            indeldata.append([ilen, l[i], iseq, index])
    return indeldata
#end

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
#end

sample_reads = {}
for sample in SAMPLES:
    sample_reads[sample] = set([i[0] for i in goodmapread(CWD + sample + '/all_good_candidate_indels.qualfilter.baq.bam')])



nums = list(map(str, list(range(10))))
canind = deque()
readl = [128, 150]
rl = deque()
rpos = deque()
for sample in SAMPLES:
    # if sample != 'WT_1': continue
    print(sample)
    smut = gooddf.loc[gooddf[7] == sample].copy()
    smut['pileuppos'] = smut[2]
    smut.loc[smut[4].str.contains('-'), 'pileuppos'] -= 1

    inds = dict(zip(zip(smut[0], smut['pileuppos']), smut[4]))
    with open(CWD + sample + '/all_good_candidate_indels.qualfilter.baq.mpileup', 'r') as fin:
        cnt = 0
        for line in tqdm(fin):
            line = line.strip().split()
            # if line[1] == '3439997':
            #     break
            if (line[0], int(line[1])) in inds:
                pinds = getindels(line[4])
                t = inds[(line[0], int(line[1]))][0]
                s = inds[(line[0], int(line[1]))][1:]
                # Positions having indels
                bind = [i for i in range(len(pinds)) if pinds[i][2].upper() != s or pinds[i][1] != t]
                for b in bind[::-1]:
                    pinds.__delitem__(b)
                # Filter if there are less than 5 reads
                if len(pinds) < 5: continue
                # Filter if >=50% reads with mutation are in the terminal 20bp
                pos = [int(line[6].split(',')[i[3]]) for i in pinds]
                flags = [line[8].split(',')[i[3]] for i in pinds]
                isr1 = [0 if format(int(i), '#014b')[7] == '1' else 1 for i in flags]
                tlen = [min(pos[i], readl[isr1[i]]-pos[i]) for i in range(len(pos))]
                rl.extend(tlen)
                rpos.extend(pos)
                # continue
                if sum([i<20 for i in tlen]) >= 0.8*len(pos): continue
                # Filter if all indel reads are from the same strand
                u = [i for i in range(len(pinds)) if pinds[i][2] == s]
                if len(u) == 0 or len(u) == len(pinds): continue
                # Filter indels which are not supported by 3 BC having at least 5 good reads
                rnames = [line[7].split(',')[p[3]] for p in pinds]
                ngrn = sum([r in sample_reads[sample] for r in set(rnames)])    # Number of good read-names
                bcs = len(set([line[9].split(',')[p[3]] for p in pinds]))
                if bcs < 3: continue
                if ngrn < 5: continue
                canind.append([sample, line[0], line[1], rnames])
with open(CWD+'canind.pickle', 'wb') as f:
    pickle.dump(canind, f)

indelreads = defaultdict(deque)
for ind in canind:
    indelreads[ind[0]].extend(ind[3])
for k in indelreads:
    indelreads[k] = set(indelreads[k])

def bamreaddata(bamfin, bamfout, pileout, reads):
    # Get bam file consisting of reads having indels
    bam = pysam.AlignmentFile(bamfin, 'rb')
    bamout = pysam.AlignmentFile(bamfout, 'wb', template=bam)
    for r in tqdm(bam):
        if r.qname in reads:
            bamout.write(r)
#end

def bamreaddata2(bamfin, bamfout, pileout, reads):
    with open(pileout, 'w') as fout:
        # for row in df.itertuples(index=False):
        run("/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools mpileup -f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta -A -q 40 -Q 0 --rf 3 --reverse-del --output-QNAME -O --output-extra RG {}".format(bamfout).split(), stdout=fout)
#end

with Pool(processes=4) as pool:
    input = [(CWD + '{s}/{s}.sorted.bt2.bam'.format(s=sample),
              CWD + '{s}/canind_reads.bt2.bam'.format(s=sample),
              CWD + '{s}/canind_reads.bt2.pileup'.format(s=sample),
              indelreads[sample]) for sample in SAMPLES]
    _ = pool.starmap(bamreaddata, input)
    _ = pool.starmap(bamreaddata2, input)

snps = defaultdict(dict)
with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out', 'r') as fin:
    for line in fin:
        line = line.strip().split()
        if line[10] != 'SNP': continue
        snps[line[0]][line[1]] = [line[3], line[4], line[9]]


sample_reads = {}
for sample in SAMPLES:
    sample_reads[sample] = {i[0]:i[1:] for i in goodmapread(CWD + sample + '/all_good_candidate_indels.qualfilter.baq.bam')}


goodind = deque()
for sample in SAMPLES:
    print(sample)
    pileout = dict()
    # Read the pileup data at SNP marker positions
    with open(CWD + '{s}/canind_reads.bt2.pileup'.format(s=sample), 'r') as fin:
        for line in tqdm(fin):
            line = line.strip().split()
            try:
                snps[line[0]][line[1]]
            except KeyError:
                continue
            else:
                pileout[(line[0], line[1])] = line[2:]

    # Get coordinates of reads
    snpreads = defaultdict(deque)
    for k, pos in pileout.items():
        for r in pos[5].split(','):
            snpreads[r].append(k)

    for ind in canind:
        if ind[0] != sample: continue
        snr = deque()
        for r in ind[3]:
            try: snpreads[r]
            except KeyError: continue
            else: snr.append(r)

        # Canind with less than three reads mapping the SNPs are filtered out
        if len(snr) < 3: continue

        # Find weather the reads has reference (R), query (Q), or NA at the SNP positions
        snrpos = deque()
        for r in snr:
            for p in snpreads[r]:
                rid = pileout[p][5].split(',').index(r)
                b = getbases(pileout[p][2])[1][rid].upper()
                if b in ['.', ',']: snrpos.append([r, p, 'R'])
                elif b == snps[p[0]][p[1]][1]: snrpos.append([r, p, 'Q'])
                else: snrpos.append([r, p, 'N'])


        if len([1 for i in snrpos if i[2] == 'R']) > 0:
            if len([1 for i in snrpos if i[2] == 'Q']) > 0: continue
            if len(set([i[0] for i in snrpos if i[2] == 'R'])) >= max(3, len(set(snr))/2):
                goodind.append(ind)
        elif len(set([i[0] for i in snrpos if i[2] == 'Q'])) >= max(3, len(set(snr))/2):
            goodind.append(ind)

selected = set([(i[0], i[1], i[2]) for i in goodind])
selind = pd.read_table(CWD + 'all_good_candidate_indels.qualfilter.bed', header=None, sep=' ')
toout = deque()
for row in selind.itertuples(index=False):
    p = row[1] if '-' in row[4] else row[2]
    if (row[7], row[0], str(p)) in selected:
        toout.append(True)
    else:
        toout.append(False)
selind['toout'] = toout
selind = selind.loc[selind['toout'] == True]
selind.drop('toout', axis=1, inplace=True)
selind.to_csv(CWD+'all_good_candidate_indels.selected.bed', index=False, header=False, sep=' ')
