# Script to analyse the identified somatic mutations
####################################################
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/software/hometools/')
from myUsefulFunctions import Namespace, extractSeq, sampfa
from collections import deque, Counter, defaultdict
import numpy as np
from scipy.stats import ttest_ind
from hometools.hometools import readfasta_iterator, writefasta


################################################################################

# <editor-fold desc="(OUTDATED) Select neighboring regions of somatic mutation and check whether they are enriched for transcription-finding motifs compared to regions without somatic mutations">

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

# </editor-fold>

################################################################################

# <editor-fold desc="(OUTDATED) Get distribution of somatic mutations in genomic regions (gene vs repeat vs TE etc)">

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

# </editor-fold>

################################################################################
# <editor-fold desc="(OUTDATED) Get distribution of somatic mutations in structural region between haplotypes (syntenic vs SR vs not_aligned/HDR) also, check overlap with genome mappability">

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
# </editor-fold>


# <editor-fold desc="Run snpEff to see the effect of mutations">

def df2vcf(df, f, ref):
    from datetime import date
    from hometools.hometools import readfasta
    with open(f, 'w') as fout:
        fout.write('##fileformat=VCFv4.3\n')
        fout.write('##fileDate=' + str(date.today()).replace('-', '') + '\n')
        fout.write('##source=syri\n')
        fout.write('\t'.join(['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']) + '\n')
        ref = readfasta(ref)
        for row in df.itertuples(index=False):
            if '-' in row[2]:
                pos = row[1]-1
                r = ref[row[0]][row[1] - 2] + row[2][1:]
                a = ref[row[0]][row[1] - 2]
            elif '+' in row[2]:
                pos = row[1]
                r = ref[row[0]][pos-1]
                a = r + row[2][1:]
            else:
                pos = row[1]
                r = ref[row[0]][pos-1]
                a = row[2]
            fout.write('\t'.join(list(map(str, [row[0], pos, '.', r, a, '.', 'PASS', '.'])))+"\n")
    return
# END


cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/sm_analysis/snpeff_run/'
refcur = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
f = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_in_all_samples.manually_selected.cleaned.vcf'
allsmmat = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_in_all_samples.manually_selected.cleaned.csv')
allsmmat = allsmmat['chromosome position alt_allele'.split()].drop_duplicates()

df2vcf(allsmmat, f, refcur)

## Run snpeff with command:
## java -jar /srv/netscratch/dep_mercier/grp_schneeberger/bin/snpEff/snpEff5.1d/snpEff.jar currot.v1 /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_sm_in_all_samples.manually_selected.cleaned.vcf > all_sm_in_all_samples.manually_selected.cleaned.ann.vcf
grep -P "LOW|MODERATE|HIGH" all_sm_in_all_samples.manually_selected.cleaned.ann.vcf  | cut -f 1,2,4,5 > affected_snps.txt

# TODO: There are 7 SMs in exons: Catalog them.
genelist = pd.read_table("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/sm_analysis/snpeff_run/snpEff_genes.txt", skiprows=1)
genelist = genelist.loc[(genelist.variants_impact_HIGH == 1) |
                        (genelist.variants_impact_LOW == 1) |
                        (genelist.variants_impact_MODERATE == 1)]

genelist = genelist[['GeneId'] + 'variants_impact_HIGH  variants_impact_LOW  variants_impact_MODERATE'.split()]

# GeneId  variants_impact_HIGH  variants_impact_LOW  variants_impact_MODERATE   variants_impact_MODIFIER
# Gene.10343                     0                    0                         1   Encodes a microtubule-associated kinase-like protein RUNKEL (RUK). Contains a putative serine/threonine kinase domain and a microtubule-binding domain. RUK directly binds to microtubules in vitro and colocalizes with mitotic preprophase band, spindle, and phragmoplast in vivo. Required for cell plate expansion in cytokinesis.
# Gene.23204                     0                    0                         1   transmembrane protein C9orf5 protein;(source:Araport11)
# Gene.26196                     0                    0                         1   NO ortholog in Thaliana. Blasts only to Prunus Dulcia
# Gene.27798                     0                    0                         1   potassium transporter
# Gene.7081                     0                    1                         0    Anion channel protein family member. Involved in negative regulation of pattern triggered immunity.
# Gene.7258                     1                    0                         0    Encodes a DNA binding protein that promotes re-arrangements of mitochondrial genome. Mutations affects mitochondrial gene expression, and impairs mitochondrial function. Dual targeting of the protein to mitochondria and chloroplasts caused by alternative translation initiation. Plastid MSH1 depletion results in variegation, abiotic stress tolerance, variable growth rate, and delayed maturity.

geneids = c('Gene.10343', 'Gene.23204', 'Gene.27798', 'Gene.7081')

# Manually checked that all these positions over-lapped exons. Further, CUR4G:6760960 is in the first two bases of a new exon.
# CUR1G	3268362	T	wt_7	L1	fruit
# CUR1G	10916855	A	mut_11_1	L1	fruit
# CUR1G	10916855	A	mut_11_2	L1	fruit
# CUR1G	11957090	G	mut_15	leaf	both
# CUR1G	11957090	G	mut_15	L2	both
# CUR3G	25500813	G	wt_19	L1	fruit
# CUR4G	6760960	G	wt_1	L1	fruit
# CUR4G	6760960	G	wt_7	L1	fruit
# CUR4G	6760960	G	wt_18	L1	fruit
# CUR4G	6760960	G	wt_19	L1	fruit
# CUR4G	6760960	G	mut_11_1	L1	fruit
# CUR4G	6760960	G	mut_11_2	L1	fruit
# CUR4G	6760960	G	mut_15	L1	fruit
# CUR4G	22366515	A	wt_7	leaf	both
# CUR4G	22366515	A	wt_7	L1	both
# CUR4G	22366515	A	wt_7	L2	both


# CUR1G	3268362	.	C	T	.	PASS	ANN=T|stop_gained|HIGH|mRNA.7258.1|Gene.7258|transcript|mRNA.7258.1|protein_coding|7/22|c.619C>T|p.Gln207*|619/3933|619/3432|207/1143||;LOF=(mRNA.7258.1|Gene.7258|1|1.00);NMD=(mRNA.7258.1|Gene.7258|1|1.00)
# CUR1G	10916855	.	C	A	.	PASS	ANN=A|missense_variant|MODERATE|mRNA.10343.1|Gene.10343|transcript|mRNA.10343.1|protein_coding|8/11|c.1030G>T|p.Val344Phe|1143/4461|1030/4101|344/1366||,A|missense_variant|MODERATE|mRNA.10343.1|Gene.10343|transcript|mRNA.10343.2|protein_coding|8/10|c.1030G>T|p.Val344Phe|1143/4263|1030/4101|344/1366||
# CUR1G	11957090	.	A	G	.	PASS	ANN=G|synonymous_variant|LOW|mRNA.7081.1|Gene.7081|transcript|mRNA.7081.1|protein_coding|22/22|c.2298A>G|p.Val766Val|2528/3560|2298/2343|766/780||,G|upstream_gene_variant|MODIFIER|mRNA.7833.1|Gene.7833|transcript|mRNA.7833.1|protein_coding||c.-1381A>G|||||1381|,G|downstream_gene_variant|MODIFIER|mRNA.9955.1|Gene.9955|transcript|mRNA.9955.1|protein_coding||c.*3218T>C|||||2994|
# CUR3G	25500813	.	A	G	.	PASS	ANN=G|missense_variant|MODERATE|mRNA.23204.1|Gene.23204|transcript|mRNA.23204.1|protein_coding|1/2|c.613A>G|p.Ile205Val|818/2520|613/2169|205/722||,G|upstream_gene_variant|MODIFIER|mRNA.23590.1|Gene.23590|transcript|mRNA.23590.1|protein_coding||c.-1880T>C|||||1880|,G|downstream_gene_variant|MODIFIER|mRNA.19298.1|Gene.19298|transcript|mRNA.19298.1|protein_coding||c.*3777A>G|||||3139|
# CUR4G	6760960	.	A	G	.	PASS	ANN=G|missense_variant&splice_region_variant|MODERATE|mRNA.27798.1|Gene.27798|transcript|mRNA.27798.1|protein_coding|8/8|c.1321A>G|p.Thr441Ala|1685/3235|1321/2370|441/789||,G|upstream_gene_variant|MODIFIER|mRNA.27664.2|Gene.27664|transcript|mRNA.27664.4|protein_coding||c.-3342A>G|||||2842|,G|upstream_gene_variant|MODIFIER|mRNA.27664.2|Gene.27664|transcript|mRNA.27664.1|protein_coding||c.-4680A>G|||||2926|,G|upstream_gene_variant|MODIFIER|mRNA.27664.2|Gene.27664|transcript|mRNA.27664.3|protein_coding||c.-4927A>G|||||2926|,G|intron_variant|MODIFIER|mRNA.27664.2|Gene.27664|transcript|mRNA.27664.2|protein_coding|2/7|c.-980-4299A>G||||||,G|intron_variant|MODIFIER|mRNA.27664.2|Gene.27664|transcript|mRNA.27664.5|protein_coding|7/8|c.657+427A>G||||||
# CUR4G	22366515	.	G	A	.	PASS	ANN=A|missense_variant|MODERATE|mRNA.26196.1|Gene.26196|transcript|mRNA.26196.1|protein_coding|2/2|c.217C>T|p.Arg73Trp|217/810|217/309|73/102||

# </editor-fold>

def getsmeffect():
    """
        THIS DID NOT RESULT IN ANY CLEAR DIFFERENCE BETWEEN WT AND SM PROTEINS

        SM affected mRNA for six transcripts. Here, I get AF2 structure for the
        mutated transcripts and then compare the structure of WT and MUT transcripts
    """
    # <editor-fold desc="Define import">
    import shutil as sh
    import pandas as pd
    import igraph as ig
    from collections import OrderedDict, defaultdict
    import numpy as np
    from sknetwork.clustering import Louvain
    from sknetwork.hierarchy import LouvainHierarchy
    from sknetwork.data import from_edge_list #, from_adjacency_list, from_graphml, from_csv
    # </editor-fold>

    # <editor-fold desc="Define defaults and constants">
    cwd = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/structurome/smeffect/'
    allpdb = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/protein_structure/currot_pdbs/all_pdbs/'
    smpdb = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/protein_structure/currot_pdbs/sm_affected_pdb/'
    snpeffdir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/sm_analysis/snpeff_run/'
    modifications = {'mRNA.7258.1':  (207, 'Q', '*'),        # Based on the list above: Gln == Q
                     'mRNA.27798.1': (441, 'T', 'A'),
                     'mRNA.10343.1': (344, 'V', 'F'),
                     'mRNA.10343.2': (344, 'V', 'F'),
                     'mRNA.23204.1': (205, 'I', 'V'),
                     'mRNA.26196.1': (73, 'R', 'W')}
    mrnaids = set(modifications.keys())

    # </editor-fold>


    # <editor-fold desc="Get sequence and structure of proteins affected by SMs">
    fin = "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.pasa_out.prot.fasta"
    affectedprot = {}
    for k, v in readfasta_iterator(open(fin, 'r'), False):
        if k in modifications:
            p = modifications[k][0] -1
            if v[p] == modifications[k][1]:
                if modifications[k][2] == '*':
                    affectedprot[k] = v[:p]
                else:
                    affectedprot[k] = v[:p] + modifications[k][2] + v[p+1:]
            else:
                print(f'ERROR: {v[p]} {modifications[k][1]}')
    writefasta(affectedprot, f'{snpeffdir}/affected_proteins.fasta')
    # Next I run alphafold to predict the structure of these mutated proteins
    # Predicted structures of mutated proteins saved here: /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/protein_structure/currot_pdbs/sm_affected_pdb/

    # </editor-fold>


    # <editor-fold desc="Get structural of proteins affected by SMs and compare them to the unmutated structure">
    # Move WT proteins and MUT proteins to working directory
    for mrna in mrnaids:
        for i in range(5):
            sh.copy(f'{allpdb}/{mrna}/ranked_{i}.pdb', f'{cwd}/pdbs/wt_{mrna}_ranked_{i}.pdb')
            sh.copy(f'{smpdb}/{mrna}/ranked_{i}.pdb', f'{cwd}/pdbs/sm_{mrna}_ranked_{i}.pdb')

    # Align the unmutated and mutated protein structures and compare them
    # Done here: SH/somatic_mutation_analysis.sh:54

    # Read the alignment and draw clusters
    aln = pd.read_table(f'{cwd}/search.m8', header=None)
    nodes = [f'{s}_{mrna}_ranked_{i}.pdb' for s in 'wt sm'.split() for mrna in mrnaids for i in range(5)]
    edges = OrderedDict()
    selectmrna = '10343'
    for row in aln.itertuples(index=False):
        if row[0] == row[1]:
            continue
        if selectmrna in row[0] and selectmrna in row[1]:
            edges[tuple(row[:2])] = 101 - row[11]


    ## Get distances for hierarchical clustering
    AM = deque()    # Adjacency matrix
    for i, s in enumerate(nodes):
        for j, e in enumerate(nodes):
            if i >= j:
                continue
            if (s, e) in edges:
                AM.append(edges[s,e])
            else:
                AM.append(1000)
    AM = list(AM)

    fig, ax = plt.subplots()
    Z = linkage(AM, method='ward')
    dendrogram(optimal_leaf_ordering(Z, AM), ax=ax, leaf_font_size=15)
    ax.set_xticks(labels=[nodes[i] for i in leaves_list(optimal_leaf_ordering(Z, AM))], ticks=plt.xticks()[0], rotation=90)
    # ax = plt.gca()
    ax.spines[:].set_visible(False)
    ax.tick_params(left=False, labelleft=False)
    # ax.set_title("Clustering of branches based on SMs")
    plt.tight_layout(pad=0.1)


    selectednodes = [n for n in nodes if selectmrna in n]
    G = ig.Graph(n=len(selectednodes))
    G.add_edges([(selectednodes.index(k[0]), selectednodes.index(k[1])) for k in edges.keys()])
    G.es['weights'] = [v for v in edges.values()]
    G.vs['sample'] = ['red' if s[:2] == 'wt' else 'blue' for s in selectednodes]
    fig, ax = plt.subplots()
    ig.plot(G, layout=G.layout_kamada_kawai(), target=ax, vertex_color=G.vs["sample"], vertex_size=1)

    adjacency = from_edge_list([(selectednodes.index(k[0]), selectednodes.index(k[1]), v)  for k,v in edges.items()])
    louvain = Louvain()
    louvain = LouvainHierarchy()
    labels = louvain.fit_predict(adjacency)
    # </editor-fold>

    return
# END
