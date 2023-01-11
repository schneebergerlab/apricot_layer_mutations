# select barcodes with high number of UMIs and sufficient numbers of reads mapping the genome
SAMPLES = ('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15')
CWD = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/barcodes/'

import pysam
from collections import deque
from multiprocessing import Pool
from functools import partial
from tqdm import tqdm
import os
import pickle
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from scipy.stats import gaussian_kde

def readbam(b, wd):
    # print(b)
    bc = b.split('-')[0]
    pysam.set_verbosity(0)
    bam = pysam.AlignmentFile(wd + b, 'rb')
    rc = 0
    um = deque()
    m = 0
    for r in bam:
        rc += 1
        try:
            um.append(r.get_tag('UB'))
        except KeyError:
            rc -= 1
            continue
        if not r.is_unmapped:
            m+=1
    return(bc, rc, len(set(um)), 0 if rc ==0 else m/rc)
#END

def getcellplots(pltname='scrna_bc_stats.pdf'):
    sampleout = {}
    with PdfPages(CWD+pltname) as pdf:
        for sample in SAMPLES:
            bams = [i for i in os.listdir(CWD + sample +'/') if '.bam' in i and sample not in i]
            with Pool(processes=62) as pool:
                bout = pool.map(partial(readbam, wd=CWD + sample +'/'), tqdm(bams))
            sampleout[sample] = bout
            # with open(CWD + sample + '/'+'bout.pickle', 'rb') as f:
            #     bout = pickle.load(f)

            bcrc = {i[0]: i[1] for i in bout}
            bcunic = {i[0]: i[2] for i in bout}
            bcmr = {i[0]: i[3] for i in bout}
            bcs = list(bcrc.keys())

            fig = plt.figure(figsize=[12, 8])
            ax = fig.add_subplot(4, 2, 1)
            ax.hist([i if i < 10000 else 10000 for i in bcrc.values()], bins=100, alpha=0.7)
            ax.set_ylabel('#BCs')
            ax.set_xlabel('Number of reads')
            ax.set_ylabel('#BCs')
            ax.minorticks_on()
            ax.grid(which='both', axis='both')
            ax.set_axisbelow(True)

            ax = fig.add_subplot(4, 2, 3)
            ax.hist(bcrc.values(), bins=100, histtype='step', cumulative=-1, linewidth=1, color='black')
            ax.minorticks_on()
            # ax.set_ylim([1, len(bcs)])
            ax.set_yscale('log')
            ax.set_xscale('log')
            ax.set_ylabel('Cumulative #BCs')
            ax.set_xlabel('Number of reads')
            ax.grid(which='both', axis='both')
            ax.set_axisbelow(True)

            ax2 = fig.add_subplot(4, 2, 2)
            ax2.hist(bcunic.values(), bins=500, histtype='step', color='black', linewidth=1)
            ax2.set_ylabel('#BCs')
            ax2.set_xlabel('Number of UMIs')
            ax2.minorticks_on()
            ax2.grid(which='both', axis='both')
            ax2.set_axisbelow(True)

            # fig=plt.figure()
            ax2 = fig.add_subplot(4, 2, 4)
            ax2.hist(bcunic.values(), bins=100, histtype='step', cumulative=-1, linewidth=1, color='black')
            # ax2.set_ylim([1, len(bcs)])
            ax2.minorticks_on()
            ax2.set_yscale('log')
            ax2.set_xscale('log')
            ax2.set_ylabel('Cumulative #BCs')
            ax2.set_xlabel('Number of UMIs')

            ax2.grid(which='both', axis='both')
            ax2.set_axisbelow(True)

            ax3 = fig.add_subplot(4, 2, 5)
            ax3.hist(bcmr.values(), bins=100, histtype='step', cumulative=True, linewidth=1, color='black')
            ax3.set_xlabel('Mapping ratio')
            ax3.set_ylabel('Cumulative #BCs')
            ax3.minorticks_on()
            ax3.grid(which='both', axis='both')
            ax3 = ax3.twinx()
            ax3.hist(bcmr.values(), bins=100, alpha=0.5)
            ax3.set_ylabel('#BCs')

            # fig = plt.figure()
            ax4 = fig.add_subplot(4, 2, 6)
            ax4.scatter([bcrc[k] for k in bcs], [bcunic[k] for k in bcs], s=0.2)
            ax4.set_xlim([1, max(bcrc.values())])
            ax4.set_ylim([1, max(bcunic.values())])
            ax4.set_xscale('log')
            ax4.set_yscale('log')
            ax4.set_xlabel('Number of reads')
            ax4.set_ylabel('Number of UMIs')

            # fig = plt.figure()
            ax4 = fig.add_subplot(4, 2, 7)
            ax4.scatter([bcrc[k] for k in bcs], [bcmr[k] for k in bcs], s=0.2)
            ax4.set_xlim([1, max(bcrc.values())])
            ax4.set_xscale('log')
            ax4.set_xlabel('Number of reads')
            ax4.set_ylabel('Mapping ratio')

            ax4 = fig.add_subplot(4, 2, 8)
            ax4.scatter([bcunic[k] for k in bcs], [bcmr[k] for k in bcs], s=0.2)
            ax4.set_xlim([1, max(bcunic.values())])
            ax4.set_xscale('log')
            # ax4.set_yscale('log')
            ax4.set_xlabel('Number of UMIs')
            ax4.set_ylabel('Mapping ratio')

            plt.suptitle(sample)
            plt.tight_layout()
            pdf.savefig(fig)
            plt.close()
    return sampleout
#END

# Get plots for
sampleout = getcellplots()

with open(CWD+'sampleout.pickle', 'rb') as f:
    sampleout = pickle.load(f)
# Number of BCs per sample:
# WT_1: 135645; WT_19: 90284; MUT_11_1: 123613; MUT_15: 149364

bcslist = {}
for sample in SAMPLES:
    bcslist[sample] = set([b[0] for b in sampleout[sample]])

# Remove BCs with UMI count <500
bcslist2 = {}
for sample in SAMPLES:
    badbcs = set([b[0] for b in sampleout[sample] if b[2] < 500])
    bcslist2[sample] = bcslist[sample] - badbcs

# Remove BCs with mapping ratio < 0.8
bcslist3 = {}
for sample in SAMPLES:
    badbcs = set([b[0] for b in sampleout[sample] if b[3] < 0.8])
    bcslist3[sample] = bcslist2[sample] - badbcs

for sample in SAMPLES:
    data = [b[2] for b in sampleout[sample] if b[0] in bcslist3[sample]]
    # data = [b[2] for b in sampleout[sample] if b[2] > 200]
    density = gaussian_kde(data)
    xs = np.linspace(0, 16000, 10000)
    density.covariance_factor = lambda: .2
    density._compute_covariance()
    plt.plot(xs, density(xs), label=sample)
plt.minorticks_on()
plt.grid(which='both', axis='both')
plt.xlabel("Number of UMI")
plt.ylabel("Frequency")
plt.legend()

# Select SAMPLE specific cutoff for minimum UMI
scut = {'WT_1': 1300, 'WT_19': 1600, 'MUT_11_1': 1200, 'MUT_15': 800}
bcslist4 = {}
for sample in SAMPLES:
    badbcs = set([b[0] for b in sampleout[sample] if b[2] < scut[sample]])
    bcslist4[sample] = bcslist3[sample] - badbcs

# Number of BCs selected:
{'WT_1': 2566, 'WT_19': 2813, 'MUT_11_1': 2965, 'MUT_15': 3047}

for sample in SAMPLES:
    with open(CWD + sample + "/" + 'good_bcs_list.txt', 'w') as fout:
        fout.write("\n".join(bcslist4[sample]))