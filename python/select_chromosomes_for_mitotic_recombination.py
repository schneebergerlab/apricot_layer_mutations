#!/usr/bin/env python3
import os
import glob
import pandas as pd
from collections import defaultdict, deque, Counter
from numpy import mean, var
from multiprocessing import Pool
from tqdm import tqdm
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.backends.backend_pdf import PdfPages
from functools import partial


CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/'
os.chdir(CWD)
SAMPLES=('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15')
INDIR='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/'
STEP = 100000
WINSIZE = 1000000

def getwindows(CHRLEN, SNPSFIN, STEP=100000, WINSIZE=1000000):
    chrs = list(CHRLEN.keys())
    snps = pd.read_table(SNPSFIN, header=None, sep='\t')
    snps = snps.drop([3, 5, 6, 7], axis=1)
    snps.columns = ['chr', 'pos', 'ref', 'alt']
    snpcnt = Counter(snps['chr'])
    binpos = defaultdict(dict)
    binlen = defaultdict(dict)
    window = defaultdict(dict)
    for c in chrs:
        cdf = snps.loc[snps['chr'] == c]
        for i in range(1, CHRLEN[c], STEP):
            binpos[c][i//STEP] = tuple(cdf.loc[(cdf['pos'] >= i) & (cdf['pos'] < i+STEP), 'pos'])
            binlen[c][i//STEP] = len(binpos[c][i//STEP])

        bins = set(binlen[c].keys())
        for i in range(1, CHRLEN[c], WINSIZE):
            binin = tuple([j//STEP for j in range(i, i+WINSIZE, STEP) if j//STEP in bins])
            nsnp = sum([binlen[c][j] for j in binin])
            if nsnp > 100: window[c][i+(WINSIZE//2)] = (binin, nsnp)
    return snps, snpcnt, binpos, binlen, window


def getcoverage_SNP(f, chrs, binpos, window):
    try:
        fdf = pd.read_table(f, header=None)
    except pd.errors.EmptyDataError:
        return
    chrvar = dict()
    for c in chrs:
        cdf = fdf.loc[fdf[0]==c]
        snprc = {row[1]: row[3] + row[5] for row in cdf.itertuples(index=False)}
        binsum = dict()
        scnt = 0    # Number of sequenced markers
        for k, pos in binpos[c].items():
            cnt = 0
            for p in pos:
                try:
                    cnt += snprc[p]
                    scnt += 1
                except KeyError: cnt += 0
            binsum[k] = cnt
        winave = dict()
        for k, bin in window[c].items():
            # winsum = sum([binsum[p] for p in bin[0]])
            try:
                winave[k] = sum([binsum[p] for p in bin[0]])/bin[1]
            except ZeroDivisionError:
                winave[k] = 0
        chrvar[c] = {'scnt': scnt,
                     'm': mean(list(winave.values())),
                     'v': var(list(winave.values()))}
    return chrvar


def readBCdata(fins, chrs, snpcnt, binpos, window):
    with Pool(processes=4) as pool:
        chrvar = pool.map(partial(getcoverage_SNP, chrs=chrs, binpos=binpos, window=window), tqdm(fins))

    chrdata = defaultdict(dict)
    for i in range(len(chrvar)):
        _, sample, bc, _ = fins[i].rsplit('/', 3)
        if chrvar[i] is not None: chrdata[sample][bc] = chrvar[i]

    df = deque()
    for sample in SAMPLES:
        for k, v in chrdata[sample].items():
            for chr in chrs:
                df.append([sample, k, chr, v[chr]['scnt'], v[chr]['m'], v[chr]['v']])
    df = pd.DataFrame(df)
    df.columns = ['sample', 'bc', 'chr', 'scnt', 'means', 'vars']
    df['vcov'] = df['vars']/df['means']
    df['nscnt'] = [row[3]/snpcnt[row[2]] for row in df.itertuples(index=False)]
    return df


def plotfigure(df, pdf, glen, title):
    samcol = {'WT_1': 'red', 'WT_19': 'blue', 'MUT_11_1': 'black', 'MUT_15': 'orange'}
    chrmar = dict(zip(list(glen), ['v', '^', '<', '>', '1', '2', '3', '4']))
    fig = plt.figure(figsize=[10, 8])
    ax1 = plt.subplot2grid((5, 5), (1, 0), rowspan=4, colspan=4)
    for chr in list(glen):
        ax1.scatter(df.loc[df['chr'] == chr, 'scnt'], df.loc[df['chr'] == chr, 'vcov'], c=df.loc[df['chr'] == chr, 'sample'].map(samcol), marker=chrmar[chr], label=chr)
    fleg = ax1.legend(ncol=2)
    plt.gca().add_artist(fleg)
    colleg = deque()
    for sample in SAMPLES:
        colleg.append(mpatches.Patch(color=samcol[sample], label=sample))
    ax1.legend(handles=colleg, loc='upper center')
    ax1.set_xlabel('Number of sequenced markers')
    ax1.set_ylabel('Variance/Mean for coverage')
    ax2 = plt.subplot2grid((5, 5), (0, 0), rowspan=1, colspan=4)
    ax2.hist(df['scnt'], bins=range(0, max(df['scnt'])+500, 500))
    ax2.set_xlabel('Number of sequenced markers')
    ax2.set_ylabel('Frequency')
    ax2.minorticks_on()
    ax2.grid(b=True, which='both', axis='x', linestyle='--', linewidth=0.25)
    ax2.set_title(title)
    ax3 = plt.subplot2grid((5, 5), (1, 4), rowspan=4, colspan=1)
    ax3.hist(df['vcov'], bins=[i/1000 for i in range(0, int((max(df['vcov']))*1000))], orientation='horizontal')
    ax3.set_ylabel('Variance/Mean for coverage')
    ax3.set_xlabel('Frequency')
    ax3.minorticks_on()
    ax3.grid(b=True, which='both', axis='y', linestyle='--', linewidth=0.25)
    plt.tight_layout()
    # plt.savefig(CWD+'coverage_variance_chromosomes.pdf')
    pdf.savefig()
    plt.close(fig)

    fig = plt.figure(figsize=[10, 8])
    ax1 = plt.subplot2grid((5, 5), (1, 0), rowspan=4, colspan=4)
    for chr in list(glen):
        ax1.scatter(df.loc[df['chr'] == chr, 'nscnt'], df.loc[df['chr'] == chr, 'vcov'], c=df.loc[df['chr'] == chr, 'sample'].map(samcol), marker=chrmar[chr], label=chr)
    ax1.set_xlim([-0.01, 0.7])
    fleg = ax1.legend(ncol=2)
    plt.gca().add_artist(fleg)
    colleg = deque()
    for sample in SAMPLES:
        colleg.append(mpatches.Patch(color=samcol[sample], label=sample))
    ax1.legend(handles=colleg, loc='upper center')
    ax1.set_xlabel('Ratio of sequenced markers')
    ax1.set_ylabel('Variance/Mean for coverage')
    ax2 = plt.subplot2grid((5, 5), (0, 0), rowspan=1, colspan=4)
    ax2.hist(df['nscnt'], bins=[i/100 for i in range(101)][:71])
    ax2.set_xlim([-0.01, 0.7])
    ax2.set_xlabel('Ratio of sequenced markers')
    ax2.set_ylabel('Frequency')
    ax2.minorticks_on()
    ax2.grid(b=True, which='both', axis='x', linestyle='--', linewidth=0.25)
    ax2.set_title(title)
    ax3 = plt.subplot2grid((5, 5), (1, 4), rowspan=4, colspan=1)
    ax3.hist(df['vcov'], bins=[i/1000 for i in range(0, int((max(df['vcov']))*1000))], orientation='horizontal')
    ax3.set_ylabel('Variance/Mean for coverage')
    ax3.set_xlabel('Frequency')
    ax3.minorticks_on()
    ax3.grid(b=True, which='both', axis='y', linestyle='--', linewidth=0.25)
    plt.tight_layout()
    pdf.savefig()
    plt.close(fig)


def writechrrc(grp, outdir):
    sample = grp[0][0]
    bc = grp[0][1]
    # Read Cur readcounts
    cfin = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/{sample}/{bc}/{sample}_{bc}_b30_q40.depth450-650.af0.4-0.6.bt2.txt'.format(sample=sample, bc=bc)
    df = pd.read_table(cfin, header=None)
    for row in grp[1].itertuples(index=False):
        cdf = df.loc[df[0] == row[2]]
        cdf.to_csv(outdir + row[2] + '/input_q40_lowvar/{chr}_{sample}_{bc}_b30_q40.depth450-650.af0.4-0.6.bt2.txt'.format(chr=row[2], sample=sample, bc=bc), sep='\t', index=False, header=False)

    ofin = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/{sample}/{bc}/{sample}_{bc}_ORA_b30_q40.depth450-650.af0.4-0.6.bt2.txt'.format(sample=sample, bc=bc)
    df = pd.read_table(ofin, header=None)
    for row in grp[1].itertuples(index=False):
        cdf = df.loc[df[0] == row[9]]
        cdf.to_csv(outdir + row[9] + '/input_q40_lowvar/{chr}_{sample}_{bc}_b30_q40.depth450-650.af0.4-0.6.bt2.txt'.format(chr=row[9], sample=sample, bc=bc), sep='\t', index=False, header=False)


CURLEN = {'CUR1G': 46975282,
          'CUR2G': 33806098,
          'CUR3G': 26861604,
          'CUR4G': 26096899,
          'CUR5G': 18585576,
          'CUR6G': 27136638,
          'CUR7G': 25539660,
          'CUR8G': 23045982}
CSNPFIN = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/strict_syn_snp_allele_readcount.depth450-650.af0.4-0.6.txt'
ORALEN = {'ORA1G': 47941709,
          'ORA2G': 31001264,
          'ORA3G': 27362024,
          'ORA4G': 28585890,
          'ORA5G': 18998592,
          'ORA6G': 27963823,
          'ORA7G': 24593044,
          'ORA8G': 22326149}
OSNPFIN = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/ora_ref_strict_syn_snp_allele_readcount.depth450-650.af0.4-0.6.txt'

CHRN = {'CUR1G': 1, 'CUR2G': 2, 'CUR3G': 3, 'CUR4G': 4, 'CUR5G': 5, 'CUR6G': 6, 'CUR7G': 7, 'CUR8G': 8, 'ORA1G': 1, 'ORA2G': 2, 'ORA3G': 3, 'ORA4G': 4, 'ORA5G': 5, 'ORA6G': 6, 'ORA7G': 7, 'ORA8G': 8}


csnps, csnpcnt, cbinpos, cbinlen, cwindow = getwindows(CURLEN, CSNPFIN, STEP=STEP, WINSIZE=WINSIZE)
osnps, osnpcnt, obinpos, obinlen, owindow = getwindows(ORALEN, OSNPFIN, STEP=STEP, WINSIZE=WINSIZE)

# Read data for CUR as reference analysis
fins = [f for sample in SAMPLES for f in glob.glob(CWD+sample+"/*/*b30_q40.depth450-650.af0.4-0.6.bt2.txt") if 'read_counts' not in f and 'ORA' not in f]
cbcdf = readBCdata(fins, list(CURLEN), csnpcnt, cbinpos, cwindow)

# Read data for ORA as reference analysis
fins = [f for sample in SAMPLES for f in glob.glob(CWD+sample+"/*/*ORA_b30_q40.depth450-650.af0.4-0.6.bt2.txt") if 'read_counts' not in f]
obcdf = readBCdata(fins, list(ORALEN), osnpcnt, obinpos, owindow)

# Save plots
with PdfPages(CWD+'coverage_variance_chromosomes.pdf') as pdf:
    plotfigure(cbcdf, pdf, CURLEN, 'Cur reference stats')
    plotfigure(obcdf, pdf, ORALEN, 'Ora reference stats')

# Filter out chromosomes with noisy sequencing based on manually selected cut-offs
maxvcov = 0.008
minscnt = 2000
minnscnt = 0.05

cbcdffilt = cbcdf.loc[(cbcdf['scnt'] >= 2000) & (cbcdf['vcov'] < 0.008) & (cbcdf['nscnt'] >= 0.05)].copy()
obcdffilt = obcdf.loc[(obcdf['scnt'] >= 2000) & (cbcdf['vcov'] < 0.008) & (obcdf['nscnt'] >= 0.05)].copy()

cbcdffilt['chrn'] = [CHRN[row[2]] for row in cbcdffilt.itertuples(index=False)]
obcdffilt['chrn'] = [CHRN[row[2]] for row in obcdffilt.itertuples(index=False)]

# Chromosomes passing cutoffs for both references
conserved = cbcdffilt.merge(obcdffilt, on=['sample', 'bc', 'chrn'], how='inner')

# Generate RTIGER input files
OUTDIR = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/rtiger_out/'
for k in CHRN:
    os.makedirs(OUTDIR + k +'/input_q40_lowvar/')
with Pool(processes=4) as pool:
    pool.map(partial(writechrrc, outdir=OUTDIR), tqdm(conserved.groupby(['sample', 'bc'])))
