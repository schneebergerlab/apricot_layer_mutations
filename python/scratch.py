import pandas as pd
from matplotlib import pyplot as plt

h1 = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_1/bam_read_counts_b30_q20.hist', header=None, sep=' ')
h2 = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_1/bam_read_counts_b30_q20.bt2.hist', header=None, sep=' ')

h1filt = h1.loc[h1[1]<=300]
h2filt = h2.loc[h2[1]<=300]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(h1filt[1],
       h1filt[0],
        color='red',
        label='minimap2: split')

ax.plot(h2filt[1],
       h2filt[0],
        color='black',
        label='bowtie2: end-to-end')
ax.set_xlabel('Coverage')
ax.set_ylabel('Frequency')
ax.legend()
fig.show()


h3 = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_1/WT_1_only_SNPs_candidate.sorted.common.bed', header=None, sep='\t')
h3filt = h3.loc[(h3[4]<=0.1) & (h3[3]>=5)]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.hist(h3filt[4],
        bins=100,
        color='blue',
        label='WT_1')

ax.set_xlabel('Variant Allele Frequency')
ax.set_ylabel('Frequency')
ax.legend()
fig.show()


h3 = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_1/WT_1_only_SNPs_candidate.sorted.bed', header=None, sep='\t')
h3filt = h3.loc[h3[4]<=0.1]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.hist(h3filt[4],
        bins=100,
        color='blue',
        label='WT_1')

ax.set_xlabel('Variant Allele Frequency')
ax.set_ylabel('Frequency')
ax.legend()
fig.show()


## Plot Strelka performance data
import pandas as pd
from matplotlib import pyplot as plt

def get_density_plot_values(data):
    import numpy as np
    from scipy.stats import gaussian_kde
    density = gaussian_kde(data)
    xs = np.linspace(min(data), max(data),1000)
    density.covariance_factor = lambda : .2
    density._compute_covariance()
    return [xs, density(xs)]


## WT_1
fn=pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/WT_1/false_negatives_by_strelka.bed', header=None, sep='\t')
fn['VAF'] = [float(i.split(';')[1].split('=')[1]) for i in fn[3]]
fn['DPR'] = [float(i.split(';')[2].split('=')[1]) for i in fn[3]]

tp=pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/WT_1/true_positive_by_strelka.bed', header=None, sep='\t')
tp['VAF'] = [float(i.split(';')[1].split('=')[1]) for i in tp[3]]
tp['DPR'] = [float(i.split(';')[2].split('=')[1]) for i in tp[3]]

fig = plt.figure(figsize=[4,5], dpi=600)
# plt.title('Performance of Strelka in calling variant')
# Plot VAF effects
ax = fig.add_subplot(2,1,1)
ax.set_xlim(0, 0.8)
fn_data = get_density_plot_values(fn['VAF'])
ax.plot(fn_data[0],
        fn_data[1],
        color='blue',
        label='False Negative')
tp_data = get_density_plot_values(tp['VAF'])
ax.plot(tp_data[0],
        tp_data[1],
        color='black',
        label='True Positive')
ax.set_xlabel('Variant Allele Frequence')
ax.set_ylabel('Frequency')
ax.legend()

# Plot DPR effect
ax = fig.add_subplot(2,1,2)
ax.set_xlim(70, 190)
fn_data = get_density_plot_values(fn['DPR'])
ax.plot(fn_data[0],
        fn_data[1],
        color='blue',
        label='False Negative')
tp_data = get_density_plot_values(tp['DPR'])
ax.plot(tp_data[0],
        tp_data[1],
        color='black',
        label='True Positive')
ax.set_xlabel('Average Depth in Region')
ax.set_ylabel('Frequency')
ax.legend()
fig.suptitle('Performance of Strelka in calling variant')
# fig.tight_layout()
plt.subplots_adjust(left=0.16, bottom=0.1, right=0.95, top=0.93, wspace=0.1, hspace=0.35)
fig.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/WT_1/strelka_performance.png')




## WT_19
fn=pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/WT_19/false_negative_by_strelka.bed', header=None, sep='\t')
fn['VAF'] = [float(i.split(';')[1].split('=')[1]) for i in fn[3]]
fn['DPR'] = [float(i.split(';')[2].split('=')[1]) for i in fn[3]]

tp=pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/WT_19/true_positive_by_strelka.bed', header=None, sep='\t')
tp['VAF'] = [float(i.split(';')[1].split('=')[1]) for i in tp[3]]
tp['DPR'] = [float(i.split(';')[2].split('=')[1]) for i in tp[3]]

fig = plt.figure(figsize=[4,5], dpi=600)
fig = plt.figure()
# plt.title('Performance of Strelka in calling variant')
# Plot VAF effects
ax = fig.add_subplot(2,1,1)
ax.set_xlim(0, 0.8)
fn_data = get_density_plot_values(fn['VAF'])
ax.plot(fn_data[0],
        fn_data[1],
        color='blue',
        label='False Negative')
tp_data = get_density_plot_values(tp['VAF'])
ax.plot(tp_data[0],
        tp_data[1],
        color='black',
        label='True Positive')
ax.set_xlabel('Variant Allele Frequence')
ax.set_ylabel('Frequency')
ax.legend()

# Plot DPR effect
ax = fig.add_subplot(2,1,2)
ax.set_xlim(70, 190)
fn_data = get_density_plot_values(fn['DPR'])
ax.plot(fn_data[0],
        fn_data[1],
        color='blue',
        label='False Negative')
tp_data = get_density_plot_values(tp['DPR'])
ax.plot(tp_data[0],
        tp_data[1],
        color='black',
        label='True Positive')
ax.set_xlabel('Average Depth in Region')
ax.set_ylabel('Frequency')
ax.legend()
fig.suptitle('Performance of Strelka in calling variant')
# fig.tight_layout()
plt.subplots_adjust(left=0.16, bottom=0.1, right=0.95, top=0.93, wspace=0.1, hspace=0.35)
fig.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/WT_1/strelka_performance.png')

top=0.935,
bottom=0.099,
left=0.099,
right=0.967,
hspace=0.255,
wspace=0.2

contigs_depth_sum = {}
with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/cur_unk.depth', 'r') as f:
    contig = ''
    sum = 0
    FIRST = True
    for line in f:
        line = line.strip().split()
        if line[0] != contig:
            if not FIRST: contigs_depth_sum[contig] = sum
            else: FIRST = False
            contig  = line[0]
            sum     = int(line[2])
            print(contig)
        else:
            sum += int(line[2])
    contigs_depth_sum[contig] = sum

df = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/cur_unk.p_ctg.fasta.chrsize', header=None)
contigs_sizes = dict(zip(df[0], df[1]))
contigs_list = list(contigs_sizes.keys())

contigs_mean_depth = {k:[contigs_depth_sum[k]/contigs_sizes[k] if k in contigs_depth_sum else 0][0] for k in contigs_list}
df['md'] = [contigs_mean_depth[k] for k in df[0]]
df.sort_values(['md'], ascending=True, inplace=True)
df['cs'] = np.cumsum(df[1])

maxcov = 100
fig = plt.figure(figsize=[6,10])
ax1 = fig.add_subplot(3,1,1)
ax1.scatter(df.md, df[1])
ax1.set_xlim([0,maxcov])
ax1.set_xlabel("mean read depth")
ax1.set_ylabel("contig size")

ax2 = fig.add_subplot(3,1,2)
plt.hist(df.md, bins=int(maxcov/2), range=[0,maxcov])
ax2.set_xlabel("mean read depth")
ax2.set_ylabel("Number of contigs")
ax2.set_xlim([0,maxcov])

ax3 = fig.add_subplot(3,1,3)
plt.plot(df.md, df.cs)
ax3.set_xlim([0, maxcov])
ax3.set_xlabel("mean read depth")
ax3.set_ylabel("Cumulative assembly size")
ax3.minorticks_on()
ax3.grid(b=True, which='minor', axis='both', linestyle='--', linewidth=0.25)
ax3.grid(b=True, which='major', axis='both', linestyle='-', linewidth=2)
plt.subplots_adjust(left=0.1, bottom=0.05, right=0.97, top=0.98, wspace=0.1, hspace=0.35)
fig.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/cur_unk_rd_stats.pdf')



contig_less2_cov = [k for k in contigs_list if contigs_mean_depth[k] > 8 or contigs_mean_depth[k] < 30]


plt.plot(range(1,(df.shape[0] + 1)), df.cs)
plt.xscale('log')
plt.yscale('log')


import pysam

f = pysam.AlignmentFile('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/good_contigs_cur_F256_F2048.sorted.bam')

a = {al.qname: al.cigarstring for al in f}
a = {al.qname: al.cigarstring for al in f if al.reference_name=='utg000114l'}
clipped = {}
clipped_count = {}
def cgtpl(cg):
    return [(int(i.split(';')[0]), i.split(';')[1]) for i in cg.replace('M',';M,').replace('I',';I,').replace('D',';D,').replace('S',';S,').replace('H',';H,').split(',')[:-1]]

for r, cg in a.items():
    if cg is None: continue
    if 'S' in cg:
        clipped[r] = cg
        clipped_count[r] = 0
        for t in cgtpl(cg):
            if t[1] in ['S', 'H']:
                clipped_count[r] += t[0]
        continue
    if 'H' in cg:
        clipped.append(cg)
        clipped_count[r] = 0
        for t in cgtpl(cg):
            if t[1] in ['S', 'H']:
                clipped_count[r] += t[0]
        continue
plt.hist(clipped_count, bins=100)

import sys
sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/')

from myUsefulFunctions import mergeRanges


df = pd.read_table('/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/tests/TMP.paf', header=None)
ctg = df[0].unique()

ctg_ovl = {}
for c in ctg:
    ctg_df = df.loc[(df[0]==c)]
    ctg_l = pd.unique(ctg_df[1])[0]
    ctg_df = ctg_df.loc[ctg_df[5]!=c]
    ctg_df_filt = ctg_df.loc[(ctg_df[9] >= 10000) & (ctg_df[9]/ctg_df[10] >= 0.9) & (ctg_df[1]<ctg_df[6])]
    rngs = mergeRanges(np.array(ctg_df_filt[[2,3]]))
    ctg_ovl[c] = [ctg_l, sum([r[1] - r[0] + 1 for r in rngs]), rngs]

select = [k for k,v in ctg_ovl.items() if v[1]/v[0] < 0.1]

f1 = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/ora/ora.v3.orasr.sorted.read_depth.stats')
f1.sort_values(['contig'], inplace=True)
f2 = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/ora/ora.v3.orapb.sorted.read_depth.stats')
f2.sort_values(['contig'], inplace=True)

from matplotlib import cm
viridis = cm.get_cmap("RdBu", 1000)
plt.scatter(f2['mean_read_depth'], f1['mean_read_depth'], c=np.log10(f1['size']), cmap=viridis, linewidth=0.2, edgecolors='black')
plt.hlines(20, xmin=0, xmax=1000)
plt.hlines(60, xmin=0, xmax=1000)
plt.vlines(10, ymin=0, ymax=1000)
plt.xscale('log')
plt.yscale('log')
plt.xlabel("PacBio coverage")
plt.ylabel("Illumina coverage")
cbar = plt.colorbar()
cbar.set_ticks(cbar.get_ticks())
cbar.set_ticklabels(np.round(10**np.array(cbar.get_ticks()), 2))
cbar.set_label("Contig Size")
plt.savefig("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/ora/orapb_orasr_read_depth.pdf")

xx = f2['mean_read_depth']
yy = f1['mean_read_depth']

import numpy as np
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
from matplotlib.colors import from_levels_and_colors
from matplotlib.collections import LineCollection
cmap, norm = from_levels_and_colors([0.0, 0.5, 1.5], ['red', 'black'])
points = np.array([xx, yy]).T.reshape(-1, 1, 2)
segments = np.concatenate([points[:-1], points[1:]], axis=1)
lines = LineCollection(segments, cmap=cmap, norm=norm)
lines.set_array(good.astype(int))
ax.add_collection(lines)
plt.show()



f1=pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/ora/ora.v3.pilon_polished.fasta.chrsize', header=None)
f2=pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/ora/ora.v3.racon_polished.fasta.chrsize', header=None)
f1.sort_values([0], inplace=True)
f2.sort_values([0], inplace=True)

plt.scatter(f1[1], f2[1])
plt.plot([0,1], [0,1], transform=plt.transAxes)
f2[2] = f2[1] - f1[1]
f2.sort_values([2], inplace=True, ascending=False)