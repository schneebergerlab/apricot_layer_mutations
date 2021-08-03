import numpy as np
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings("error")

snps = defaultdict(dict)
with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out', 'r') as fin:
    for line in fin:
        line = line.strip().split()
        if line[10] != 'SNP': continue
        if line[3] == 'N' or line[4] == 'N': continue
        snps[line[0]][line[1]] = [line[3], line[4], line[9]]

garb = dict()
def snpreadcnt(fin, snps):
    BASEDICT={'A': 4, 'C': 5, 'G': 6, 'T': 7}
    sdata = deque()
    with open(fin, 'r') as fin:
        for line in tqdm(fin):
            line = line.strip().split()
            if int(line[3]) < 80: continue
            if int(line[3]) > 250: continue
            try: snps[line[0]][line[1]]
            except KeyError: pass
            else: sdata.append((line[0], line[1], int(line[BASEDICT[snps[line[0]][line[1]][0]]]), int(line[BASEDICT[snps[line[0]][line[1]][1]]])))
    return sdata

with Pool(processes=4) as pool:
    fins = [CWD + '{s}/bam_read_counts_b30_q10.bt2.txt'.format(s=s) for s in SAMPLES]
    garb = pool.map(partial(snpreadcnt, snps=snps), fins)

with open('TMP_snpreadcnt.pickle', 'rb') as f:
    garb = pickle.load(f)

snpmaf = defaultdict(deque)
for i in range(4):
    for s in garb[i]:
        if s[2] + s[3] == 0:
            raise ValueError(s)
        snpmaf[(s[0], s[1])].append(s[3]/(s[2]+s[3]))

# snpvar = [np.var(v)/np.mean(v) for v in snpmaf.values()]
snpoutlierdist = deque()
for k, v in tqdm(snpmaf.items()):
    if len(v) < 4: continue
    try:
        snpoutlierdist.append(min(v - np.mean(v))/np.mean(v))
    except RuntimeWarning:
        snpoutlierdist.append(0)
        # break
plt.hist(snpoutlierdist, bins=100)
plt.ylabel("Number of SNPs")
plt.xlabel('Normalised max difference from mean MAF')
plt.tight_layout()
plt.savefig(CWD+'mean_difference_histogram_homqry.png', dpi=600)



snpoutlierdist2 = deque()
for k, v in tqdm(snpmaf.items()):
    if len(v) < 4: continue
    try:
        snpoutlierdist2.append(max(v - np.mean(v))/np.mean(v))
    except RuntimeWarning:
        snpoutlierdist2.append(0)
plt.hist(snpoutlierdist2, bins=100)
plt.ylabel("Number of SNPs")
plt.xlabel('Normalised max difference from mean MAF')
plt.tight_layout()
plt.savefig(CWD+'mean_difference_histogram_homref.png', dpi=600)



for i in range(4):
    plt.scatter([j[0] for j in garb[i]], [j[1] for j in garb[i]], alpha=0.25, label=SAMPLES[i], s=0.1)


