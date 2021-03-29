from collections import defaultdict, deque

SAMPLES = set(['WT_1', 'WT_19', 'MUT_11_1', 'MUT_15'])
scan = defaultdict(deque)

can_reg = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/common_candidates.regions'
with open(can_reg, 'r') as fin:
    for line in fin:
        line = line.strip().split()
        scan[(SAMPLES - set(line[4].split(','))).pop()].append(line[0] + '_' + line[2])

mutpos = [v for values in scan.values() for v in values]

smut = defaultdict(dict)
for sample in SAMPLES:
    with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/{s}/{s}.gatk_hc.snps.bt2.candidates.bed'.format(s=sample), 'r') as fin:
        for line in fin:
            line = line.strip().split()
            smut[sample][line[0] + '_' + line[2]] = (line[5], line[6])

# validate the mutations
badpos = deque()
muts = {}
for mut in mutpos:
    c = []
    for sample in SAMPLES:
        if mut in smut[sample]: c.append(smut[sample][mut])
    if len(set(c)) != 1:
        badpos.append(mut)
    else:
        muts[mut] = list(set(c))[0]

for sample in SAMPLES:
    for pos in badpos:
        try:
            scan[sample].remove(pos)
        except:
            pass
        if pos in smut[sample]: smut[sample].pop(pos)

for pos in badpos:
    mutpos.remove(pos)

BASE_DICT = {'A': 4, 'C': 5, 'G': 6, 'T': 7}

scnt = defaultdict(dict)
for sample in SAMPLES:
    with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/{}/common_candidates.read_count.txt'.format(sample), 'r') as fin:
        for line in fin:
            line = line.strip().split()
            if line[0] + '_' + line[1] in badpos: continue
            scnt[sample][line[0] + '_' + line[1]] = [int(line[BASE_DICT[muts[line[0] + '_' + line[1]][0]]]), int(line[BASE_DICT[muts[line[0] + '_' + line[1]][1]]])]

[print(s, k ) for s, values in scnt.items() for k, v in values.items() if v[1] < 10]

mutAF = {}
for mut in muts:
    af = deque()
    for sample in SAMPLES:
        af.append(scnt[sample][mut][1]/(scnt[sample][mut][0]+scnt[sample][mut][1]))
    mutAF[mut] = af

import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot()
for mut, v in mutAF.items():
    if all([m > 0.1 for m in v]) : continue
    if v[3] < 0.3 and v[2] > 0.3:
        print(mut, v)
    ax.plot([1,2,3,4], v, label=mut)

ax.legend()






