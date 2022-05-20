# Script to analyse the identified somatic mutations
####################################################
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/software/hometools/')
from myUsefulFunctions import Namespace, extractSeq, sampfa
from collections import deque, Counter
import numpy as np
from scipy.stats import ttest_ind

def getsmpos():
    leafsm = pd.read_table("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/high_cov_mutants_sorted.all_samples.selected.txt", header=None)
    layersm = pd.read_table("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/high_conf_layer_specific_somatic_mutations.selected.txt")
    layersm.columns = [0, 1, 2, 3] + list(layersm.columns[4:])
    CWD = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/sm_analysis/smlist_20052022/'
    BUFF = 5000
    REFCUR = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
    leafpos = leafsm[[0, 1]].copy()
    layerpos = layersm[[0, 1]].copy()
    smpos = pd.concat([leafpos, layerpos])
    smpos.sort_values([0, 1], inplace=True)
    smpos.reset_index(drop=True, inplace=True)
    smpos.drop_duplicates(inplace=True)
    return smpos
#END

################################################################################
# Select neighboring regions of somatic mutation and check whether they are enriched for transcription-finding motifs compared to regions without somatic mutations

## Get somatic-mutation locations and 5kb region around them : Read the somatic mutaion location in leafs and the layer-specific somatic mutations

# Extract sequence at somatic mutation regions
with open(f"{CWD}/sm_reg_coords.txt", "w") as fout:
    for row in smpos.itertuples(index=False):
        fout.write("\t".join(list(map(str, [row[0], row[1]-int(BUFF/2)-1, row[1]+int(BUFF/2)]))) + "\n")
args = Namespace(chr=None, start=None, end=None,
                 fasta=Namespace(name=REFCUR),
                 o=Namespace(name=f"{CWD}/sm_reg.fa"),
                 fin=Namespace(name=f"{CWD}/sm_reg_coords.txt"))
extractSeq(args)
# sample regions from the chromosomes
args = Namespace(fa=Namespace(name=REFCUR),
                 n=100, s=BUFF,
                 chr=['CUR1G', 'CUR2G', 'CUR3G', 'CUR4G', 'CUR5G', 'CUR6G', 'CUR7G', 'CUR8G'],
                 o=Namespace(name=f"{CWD}/gen_samp.fa"),
                 v=None)
sampfa(args)

## Run fimo: commands in somatic_mutation_analysis.sh

## Get distribution of TFB motifs in the regions and the distribution of TFB overlapping the SM/center of region
sm_motifs = deque()
with open(f"{CWD}fimo_sm/fimo.tsv", 'r') as fin:
    line = fin.readline()
    for line in fin:
        if line[0] == '#': continue
        line = line.strip().split()
        if line == []: continue
        sm_motifs.append([line[0], line[2], line[3], line[4], line[8], line[9]])
sm_motifs = pd.DataFrame(sm_motifs)
sm_motifs[[2, 3]] = sm_motifs[[2, 3]].astype(int)

sample_motifs = deque()
with open(f"{CWD}fimo_sample/fimo.tsv", 'r') as fin:
    line = fin.readline()
    for line in fin:
        if line[0] == '#': continue
        line = line.strip().split()
        if line == []: continue
        sample_motifs.append([line[0], line[2], line[3], line[4], line[8], line[9]])
sample_motifs = pd.DataFrame(sample_motifs)
sample_motifs[[2, 3]] = sample_motifs[[2, 3]].astype(int)

## Draw summary plots
fig = plt.figure()
### Distribution of number of motifs per sequence
ax = fig.add_subplot(2, 2, 1)
sm_motif_cnt = np.array(list(Counter(sm_motifs[1]).values()))
sample_motif_cnt = np.array(list(Counter(sample_motifs[1]).values()))
ax.violinplot([sm_motif_cnt, sample_motif_cnt], positions=[1, 2])
ax.set_xticks(ticks=[1, 2], labels=["sm", "control"])
ax.set_ylabel("Total number of motifs")

### Distribution of number of unique motifs per sequence
ax = fig.add_subplot(2, 2, 2)
g = sm_motifs.drop_duplicates(subset=[0, 1])
sm_motif_cnt = np.array(list(Counter(g[1]).values()))
g = sample_motifs.drop_duplicates(subset=[0, 1])
sample_motif_cnt = np.array(list(Counter(g[1]).values()))
ax.violinplot([sm_motif_cnt, sample_motif_cnt], positions=[1, 2])
ax.set_xticks(ticks=[1, 2], labels=["sm", "control"])
ax.set_ylabel("Number of unique motifs")
ttest_ind(sm_motif_cnt, sample_motif_cnt, alternative="greater")

### Distribution of number of unique overlapping the SM/center of region
ax = fig.add_subplot(2, 2, 3)
g = sm_motifs.drop_duplicates(subset=[0, 1])
g = g.loc[(g[2] <= 2500) & (g[3] >= 2500)]
sm_motif_cnt = np.array(list(Counter(g[1]).values()))
g = sample_motifs.drop_duplicates(subset=[0, 1])
g = g.loc[(g[2] <= 2500) & (g[3] >= 2500)]
sample_motif_cnt = np.array(list(Counter(g[1]).values()))
ax.violinplot([sm_motif_cnt, sample_motif_cnt], positions=[1, 2])
ax.set_xticks(ticks=[1, 2], labels=["sm", "control"])
ax.set_ylabel("Overlapping motif count")
ttest_ind(sm_motif_cnt, sample_motif_cnt, alternative="greater")
plt.suptitle("Stats for transcription binding motifs in SM regions and randomly selected regions")
plt.tight_layout()
################################################################################

# Get distribution of somatic mutations in genomic regions (gene vs repeat vs TE etc)
smpos = getsmpos()
