# Script to analyse the identified somatic mutations
####################################################
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/software/hometools/')
from myUsefulFunctions import Namespace, extractSeq, sampfa, readfasta
from collections import deque, Counter, defaultdict
import numpy as np
from scipy.stats import ttest_ind

def getsmpos():
    leafsm = pd.read_table("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/high_cov_mutants_sorted.all_samples.selected.txt", header=None)
    layersm = pd.read_table("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/high_conf_layer_specific_somatic_mutations.selected.txt")
    layersm.columns = [0, 1, 2, 3] + list(layersm.columns[4:])
    leafpos = leafsm[[0, 1]].copy()
    layerpos = layersm[[0, 1]].copy()
    smpos = pd.concat([leafpos, layerpos])
    smpos.sort_values([0, 1], inplace=True)
    smpos.reset_index(drop=True, inplace=True)
    smpos.drop_duplicates(inplace=True)
    smpos.reset_index(drop=True, inplace=True)
    return smpos
#END

def df2vcf(df, f, ref):
    from datetime import date
    df = smpos.drop_duplicates(['chr', 'pos', 'ref', 'alt'])
    with open(f, 'w') as fout:
        fout.write('##fileformat=VCFv4.3\n')
        fout.write('##fileDate=' + str(date.today()).replace('-', '') + '\n')
        fout.write('##source=syri\n')
        fout.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']) + '\n')
        ref=readfasta(ref)
        for row in df.itertuples(index=False):
            if '-' in row[3]:
                pos = row[1]-1
                r = ref[row[0]][row[1] - 2] + row[3][1:]
                a = ref[row[0]][row[1] - 2]
            elif '+' in row[3]:
                pos = row[1]
                r = row[2]
                a = row[2] + row[3][1:]
            else:
                pos = row[1]
                r = row[2]
                a = row[3]
            fout.write('\t'.join(list(map(str, [row[0], pos, '.', r, a, '.', 'PASS', '.'])))+"\n")
#END

################################################################################
# Select neighboring regions of somatic mutation and check whether they are enriched for transcription-finding motifs compared to regions without somatic mutations

## Get somatic-mutation locations and 5kb region around them : Read the somatic mutaion location in leafs and the layer-specific somatic mutations
CWD = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/sm_analysis/smlist_20052022/'
# BUFF = 5000
BUFF = 20
REFCUR = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
smpos = getsmpos()
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
g = g.loc[(g[2] <= int(BUFF/2)) & (g[3] >= int(BUFF/2))]
sm_motif_cnt = np.array(list(Counter(g[1]).values()))
g = sample_motifs.drop_duplicates(subset=[0, 1])
g = g.loc[(g[2] <= int(BUFF/2)) & (g[3] >= int(BUFF/2))]
sample_motif_cnt = np.array(list(Counter(g[1]).values()))
ax.violinplot([sm_motif_cnt, sample_motif_cnt], positions=[1, 2])
ax.set_xticks(ticks=[1, 2], labels=["sm", "control"])
ax.set_ylabel("Overlapping motif count")
ttest_ind(sm_motif_cnt, sample_motif_cnt, alternative="greater")
plt.suptitle(f"Stats for transcription binding motifs in SM regions and randomly selected regions. Buff:{BUFF}")
plt.tight_layout()

################################################################################
# Get distribution of somatic mutations in genomic regions (gene vs repeat vs TE etc)

## Read somatic mutation list
leafsm = pd.read_table("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/high_cov_mutants_sorted.all_samples.selected.txt", header=None)
layersm = pd.read_table("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/high_conf_layer_specific_somatic_mutations.selected.txt")
layersm.columns = [0, 1, 2, 3] + list(layersm.columns[4:])
leafsm = leafsm[[0, 1, 2, 3, 6]]
leafsm.columns = ['chr', 'pos', 'ref', 'alt', 'sample']
leafsm['layer'] = '-'
layersm = layersm[[0, 1, 2, 3, 'Layer']]
layersm.columns = ['chr', 'pos', 'ref', 'alt', 'layer']
layersm['sample'] = 'MUT_11_1'
smpos = pd.concat([leafsm, layersm])
smpos.sort_values(['chr', 'pos'], inplace=True)
smpos.reset_index(drop=True, inplace=True)
smpos.reset_index(drop=True, inplace=True)


repeatlist = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/repeat/RepeatMasker/cur.genome.v1.fasta.out'
tegenes = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/EVM_PASA/pasa_on_mancur/cur.pasa_out.sort.TE.gff3'
progenes = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/EVM_PASA/pasa_on_mancur/cur.pasa_out.sort.protein_coding.3utr.gff3'

## Read repeat region coordinates
repcoord = deque()
with open(repeatlist, 'r') as fin:
    for i in range(3):
        line = fin.readline()
    for line in fin:
        line = line.strip().split()
        repcoord.append([line[4], line[5], line[6], line[8], line[10]])
repcoord = pd.DataFrame(repcoord)
repcoord[[1, 2]] = repcoord[[1, 2]].astype(int)
repcoord[5] = 'repeat'

## Read TE gene coordinates
tecoord = deque()
with open(tegenes, 'r') as fin:
    for line in fin:
        line = line.strip().split()
        if line[2] == 'gene':
            tecoord.append([line[0], line[3], line[4], line[6], line[8].split(";")[0].rsplit("=")[1]])
tecoord = pd.DataFrame(tecoord)
tecoord[[1, 2]] = tecoord[[1, 2]].astype(int)
tecoord[5] = 'tegenes'

## Read protein coding gene coordinates
procoord = deque()
with open(progenes, 'r') as fin:
    for line in fin:
        line = line.strip().split()
        if line[2] == 'gene':
            procoord.append([line[0], line[3], line[4], line[6], line[8].split(";")[0].rsplit("=")[1]])
procoord = pd.DataFrame(procoord)
procoord[[1, 2]] = procoord[[1, 2]].astype(int)
procoord[5] = 'progenes'

annocoord = pd.concat([repcoord, tecoord, procoord])
annocoord.columns = ['Chromosome', 'Start', 'End', 'Strand', 'ID', 'type']
annocoord.loc[annocoord['Strand'] == 'C', 'Strand'] = '-'

smanno = deque()
for grp in smpos.groupby(['chr']):
    chrcoord = annocoord.loc[annocoord['Chromosome'] == grp[0]]
    for row in grp[1].itertuples(index=False):
        df = chrcoord.loc[(chrcoord['Start'] <= row[1]) & (chrcoord['End'] >= row[1])].copy()
        df['chr'] = row[0]
        df['pos'] = row[1]
        smanno.append(df)
smanno = pd.concat(smanno)
smanno[['Start', 'End']] = smanno[['Start', 'End']].astype(str)
smanno.drop_duplicates(inplace=True)
smanno = smpos.merge(smanno, how='left', on=['chr', 'pos'])
smanno.fillna('-', inplace=True)
smanno.to_csv("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/somatic_mutation_genomic_region_annotation.txt", index=False, sep='\t')

################################################################################
# Get distribution of somatic mutations in structural region between haplotypes (syntenic vs SR vs not_aligned/HDR) also, check overlap with genome mappability

## Read smpos same way as in repeat/gene overlap section above
syripath = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/all_sv.syri.out'
syriout = deque()
with open(syripath, 'r') as fin:
    for line in fin:
        line = line.strip().split()
        if line[10] in {'SNP', 'INS', 'DEL'}: continue
        syriout.append(line)
syriout = pd.DataFrame(syriout)
syriout = syriout.loc[~(syriout[0] == '-')]
syriout = syriout.loc[~(syriout[10].isin(['SYN', 'INV', 'TRANS', 'DUP', 'INVTR', 'INVDP']))]
syriout[[1, 2]] = syriout[[1, 2]].astype(int)
syriout[9] = syriout[9].str.replace(r'[0-9]', '', regex=True)
smanno = deque()
for grp in smpos.groupby(['chr']):
    chrcoord = syriout.loc[syriout[0] == grp[0]]
    for row in grp[1].itertuples(index=False):
        df = chrcoord.loc[(chrcoord[1] <= row[1]) & (chrcoord[2] >= row[1])].copy()
        df['chr'] = row[0]
        df['pos'] = row[1]
        smanno.append(df)
smanno = pd.concat(smanno)
smanno = smanno[[0, 1, 2, 8, 9, 10, 11, 'chr', 'pos']]
smanno[[1, 2]] = smanno[[1, 2]].astype(str)
smanno.drop_duplicates(inplace=True)
smanno = smpos.merge(smanno, how='left', on=['chr', 'pos'])
smanno.fillna('-', inplace=True)
smanno.columns = list(smanno.columns)[:6] + ['chrom', 'start', 'end', 'id', 'parentid', 'type', 'duptype']
smanno.to_csv("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/somatic_mutation_structural_annotation.txt", index=False, sep='\t')

################################################################################
# Get correlation of mappability and SM. Ideally, this should not be different
# than the background distribution

## Read smpos same way as in repeat/gene overlap section above
genmap = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.genmap.E0_K51.map.bedgraph'
smindex = 0
smanno1 = deque()
c, p = -1, -1
with open(genmap, 'r') as fin:
    while smindex < smpos.shape[0]:
        while tuple(smpos.iloc[smindex, [0, 1]]) == (c, p):
            smindex += 1
        c, p = smpos.iloc[smindex, [0, 1]]
        for line in fin:
            line = line.strip().split()
            if line[0] != c: continue
            if int(line[1])-25 > p: continue
            if int(line[2])-25 < p: continue
            smanno1.append((c, p, round(float(line[3]), 2)))
            break
        smindex += 1
genmap = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.genmap.E1_K51.map.bedgraph'
smindex = 0
smanno2 = deque()
c, p = -1, -1
with open(genmap, 'r') as fin:
    while smindex < smpos.shape[0]:
        while tuple(smpos.iloc[smindex, [0, 1]]) == (c, p):
            smindex += 1
        c, p = smpos.iloc[smindex, [0, 1]]
        for line in fin:
            line = line.strip().split()
            if line[0] != c: continue
            if int(line[1])-25 > p: continue
            if int(line[2])-25 < p: continue
            smanno2.append((c, p, round(float(line[3]), 2)))
            break
        smindex += 1
smanno1 = pd.DataFrame(smanno1)
smanno1.columns = ['chr', 'pos', 'genmap_E0']
smanno2 = pd.DataFrame(smanno2)
smanno2.columns = ['chr', 'pos', 'genmap_E1']
smanno = smanno1.merge(smanno2, how='left', on=['chr', 'pos'])
smanno = smpos.merge(smanno, how='left', on=['chr', 'pos'])
smanno.fillna('-', inplace=True)
smanno.to_csv("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/somatic_mutation_genome_mappability.txt", index=False, sep='\t')

################################################################################
# Run snpEff to see the effect of mutations
cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/sm_analysis/smlist_20052022/snpeff_run/'
f = f'{cwd}sm.vcf'
df2vcf(smpos, f, REFCUR)

################################################################################
# Check whether the mutations are in CpG and dipyrimidine sites
ref = readfasta(REFCUR)
SNPs = smpos.loc[~smpos['alt'].str.contains('-|\+', regex=True)].copy()
SNPs.drop_duplicates(['chr', 'pos', 'alt'], inplace=True)
seq = deque()
for row in SNPs.itertuples(index=False):
    seq.append((ref[row[0]][row[1]-2: row[1]], ref[row[0]][row[1]-1: row[1]+1]))

################################################################################
################################################################################
################################################################################

## ReDo analysis with the mutations from all leaves and layer data (30.01.2023)
## Read somatic mutation list
sdict = dict(zip(('WT_1', 'wt7', 'wt18', 'WT_19', 'MUT_11_1', 'mut11_2', 'mut4', 'MUT_15'), ('wt_1', 'wt_7', 'wt_18', 'wt_19', 'mut_11_1', 'mut_11_2', 'mut_4', 'mut_15')))
refcur = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'

leafsm = pd.read_table("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/high_cov_mutants_sorted.all_samples.selected.txt", header=None)
layersm = pd.read_table("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/all_layer_somatic_variants.filtered.txt")
leafsm = leafsm[[0, 1, 2, 3, 6]]
leafsm.columns = ['chromosome', 'position', 'ref_allele', 'alt_allele', 'branch']
leafsm['Layer'] = '-'
layersm = layersm[['chromosome', 'position', 'ref_allele', 'alt_allele', 'branch', 'Layer']]
leafsm['branch'] = [sdict[b] for b in leafsm['branch']]

smpos = pd.concat([leafsm, layersm])
smpos.sort_values(['chromosome', 'position'], inplace=True)
smpos.reset_index(drop=True, inplace=True)
f ='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/sm_analysis/smlist_30012023/somatic_mutations.vcf'
df2vcf(smpos, f, refcur)

def df2vcf(df, f, ref):
    from datetime import date
    from hometools.hometools import readfasta
    df = df.drop_duplicates(['chromosome', 'position', 'ref_allele', 'alt_allele'])
    with open(f, 'w') as fout:
        fout.write('##fileformat=VCFv4.3\n')
        fout.write('##fileDate=' + str(date.today()).replace('-', '') + '\n')
        fout.write('##source=syri\n')
        fout.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']) + '\n')
        ref=readfasta(ref)
        print(df.shape)
        for row in df.itertuples(index=False):
            if '-' in row[3]:
                pos = row[1]-1
                r = ref[row[0]][row[1] - 2] + row[3][1:]
                a = ref[row[0]][row[1] - 2]
            elif '+' in row[3]:
                pos = row[1]
                r = row[2]
                a = row[2] + row[3][1:]
            else:
                pos = row[1]
                r = row[2]
                a = row[3]
            fout.write('\t'.join(list(map(str, [row[0], pos, '.', r, a, '.', 'PASS', '.'])))+"\n")
#END
## Run snpeff with command:
## java -jar /srv/netscratch/dep_mercier/grp_schneeberger/bin/snpEff/snpEff5.1d/snpEff.jar currot.v1 somatic_mutations.vcf > somatic_mutations.ann.vcf

# Get SMs in mutant branches
smmut = smpos.loc[smpos.branch.isin(['mut_11_1', 'mut_11_2', 'mut_15'])]
samplecnt = {grp[0]: grp[1].shape[0] for grp in smmut.groupby(['chromosome', 'position', 'alt_allele'])}
smmut.drop_duplicates(['chromosome', 'position', 'ref_allele', 'alt_allele'], inplace=True)
# Read snpeff output VCF and select mutant SMs
# pos = set(zip(smmut.chromosome, smmut.position))
pos = set([(k[0], k[1]) for k, v in samplecnt.items() if v == 3])
posann = deque()
with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/sm_analysis/smlist_30012023/snpeff_run/somatic_mutations.ann.vcf', 'r') as fin:
    for line in fin:
        if line[0] == '#': continue
        line = line.strip().split()
        if (line[0], int(line[1])) in pos:
            posann.append(line)
# Get list of affected/nearby genes
genesann = deque()
for p in posann:
    a = p[7].split('=')[1].split(',')
    for a1 in a:
        b = a1.split('|')[3]
        if b == '': continue
        if '-' in b:
            genesann.extend(b.split('-'))
        else:
            genesann.append(b)
genesann = set(genesann)
with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/sm_analysis/smlist_30012023/snpeff_run/somatic_mutations.affected_genes.txt', 'w') as fout:
    fout.write('\n'.join(genesann))
