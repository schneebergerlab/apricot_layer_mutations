#!usr/bin/env python
"""
Filter the COs predicted by RTIGER. Filtering is done based on size, as well as
Allele Frequency.
Finally, filtered COs (homozygous regions) predicted when using both genomes are
considered as TRUE homozygous regions.
"""

# Import libraries
import os
import re
import pandas as pd
from collections import deque, defaultdict
from matplotlib import pyplot as plt
import numpy as np
from multiprocessing import Pool
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import ttest_ind
import sys
sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/software/hometools/')
from myUsefulFunctions import mergepdf
from myUsefulFunctions import readfasta
from intervaltree import IntervalTree, Interval
import pickle
sys.path.insert(0, '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/')
from select_chromosomes_for_mitotic_recombination import getwindows
from functools import partial
from tqdm import tqdm
################################################################################
def getsynregions(f):
    synreg = deque()
    with open(f, 'r') as fin:
        for line in fin:
            if 'SYN' not in line: continue
            if 'SYNAL' in line: continue
            if 'SNP' in line: continue
            line = line.strip().split()
            if line[10] == 'SYN':
                synreg.append((line[0], int(line[1]), int(line[2]), line[5], int(line[6]), int(line[7])))
    return synreg


def getcelldata(df, sample, bc, chr):
    if sample != 'NA':
        garb = df.loc[(df['sample'] == sample)].copy()
    else:
        garb = df.copy()
    if bc != 'NA':
        garb = garb.loc[(garb['bc'] == bc)]
    if chr != 'NA':
        garb = garb.loc[(garb['chr'] == chr)]
    return garb


def readco(codir):
    print(codir)
    chrdf = deque()
    # os.chdir(CWD + chr + '/rtiger_co_q40_filter/')     # Analyse samples with q=40
    os.chdir(codir)     # Analyse samples with q=40
    dirs = [d for d in os.listdir() if os.path.isdir(d)]
    for d in dirs:
        s = [s for s in SAMPLES if s+'_' in d][0]
        # print(s)
        bc = re.findall(r'[a-zA-Z]{16}', d)[0]
        g = [f for f in os.listdir(d) if 'CompleteBlock-state' in f][0]
        bcdf = pd.read_table(d + '/' + g, header=None)
        bcdf.columns = ['chr', 'start', 'end', 'genotype']
        bcdf['sample'] = s
        bcdf['bc'] = bc
        bcdf = bcdf[['sample', 'bc', 'chr', 'start', 'end', 'genotype']]
        # df = pd.concat([df, bcdf])
        chrdf.append(bcdf)
    return pd.concat(chrdf)


def readterm(chr):
    dirpath = CWD + chr + "/input_q40_filter/"
    fins = os.listdir(dirpath)
    chrdata = deque()
    for f in fins:
        s = [s for s in SAMPLES if s+'_' in f][0]
        bc = re.findall(r'[a-zA-Z]{16}', f)[0]
        cnt = 0
        with open(dirpath + f, 'r') as fin:
            for line in fin:
                cnt += 1
                if cnt < 101: continue
                else:
                    chrdata.append([s, bc, chr, int(line.strip().split()[1])])
                    break
    return pd.DataFrame(chrdata)


def mergegenotype(indf):
    outdf = pd.DataFrame()
    for grp in indf.groupby(['sample', 'bc', 'chr', 'ref']):
        # if len(pd.unique(grp[1]['genotype'])) == 1: continue
        tmp = deque()
        current = ''
        for row in grp[1].itertuples(index=False):
            if current == '':
                current = list(row)
                continue
            if current[5] != row.genotype:
                tmp.append(current)
                current = list(row)
            else:
                current[4] = row.end
        tmp.append(current)
        tmp = pd.DataFrame(tmp)
        tmp.columns = grp[1].columns
        outdf = pd.concat([outdf, tmp])
    return outdf


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


def getcocells(grp, insnp, fins):
    """
    grp: df corresponding to a BC chromosome
    snps: snp position in the chromosome
    """
    # print(grp[0])
    if grp[1]['highhom'][0] == 1:
        return
    if len(pd.unique(grp[1]['genotype'])) == 1:
        return
    # snps = defaultdict(list)
    # for row in SNPS.loc[SNPS['chr']==grp[0][2]].itertuples(index=False):
    #     snps[row.pos] = [0, 0]
    # with open('{INDIR}/rtiger_out/{chr}/input_q40_filter/{chr}_{s}_{bc}_b30_q40.depth350-550.af0.35-0.6.bt2.txt'.format(INDIR=INDIR, s=grp[0][0], bc=grp[0][1], chr=grp[0][2]), 'r') as fin:
    snps = insnp[grp[0][2]]
    with open(fins.format(s=grp[0][0], bc=grp[0][1], chr=grp[0][2]), 'r') as fin:
        for line in fin:
            line = line.strip().split()
            snps[int(line[1])] = [int(line[3]), int(line[5])]
    snppos = np.array(list(snps.keys()))
    hetpos = deque([])
    for row in grp[1].itertuples(index=False):
        if row.genotype in ['AA', 'BB']: continue
        hetpos.extend(snppos[(snppos >= row.start) & (snppos < row.end)])
    if len(hetpos) < 1000:
        print("Too little het markers {}".format(grp[0]))
        return
    hetpos = np.array(hetpos)
    rcnthet = np.array([snps[h][0] for h in hetpos])
    qcnthet = np.array([snps[h][1] for h in hetpos])
    rmeanhet = [np.mean(np.random.choice(rcnthet, 500, replace=False)) for _ in range(100)]
    qmeanhet = [np.mean(np.random.choice(qcnthet, 500, replace=False)) for _ in range(100)]

    cnt = 0
    newgen = deque()

    for row in grp[1].itertuples(index=False):
        # cnt += 1
        # if cnt == 1: break
        if row.genotype == 'AB':
            newgen.append('AB')
            continue
        hompos = snppos[(snppos >= row.start) & (snppos < row.end)]
        if row.genotype == 'AA':
            h = rmeanhet
            l = qmeanhet
            hcnthom = np.array([snps[h][0] for h in hompos])
            lcnthom = np.array([snps[h][1] for h in hompos])
        else:
            h = qmeanhet
            l = rmeanhet
            hcnthom = np.array([snps[h][1] for h in hompos])
            lcnthom = np.array([snps[h][0] for h in hompos])
        if sum(hcnthom>0) < 250:
            # print(Counter(hcnthom))
            newgen.append('AB')
            continue
        hmeanhom = [np.mean(np.random.choice(hcnthom, 500, replace=False)) for _ in range(100)]
        lmeanhom = [np.mean(np.random.choice(lcnthom, 500, replace=False)) for _ in range(100)]
        if np.mean(lmeanhom) > 0.5*np.mean(hmeanhom):
            newgen.append('AB')
        elif ttest_ind(hmeanhom, h, alternative='greater', equal_var=False).pvalue > 0.001:
            newgen.append('AB')
        # elif np.mean(hmeanhom) < 1.5*np.mean(h) or np.mean(hmeanhom) > 3*np.mean(h):
        elif np.mean(hmeanhom) < 1.5*np.mean(h):
            newgen.append('AB')
        elif ttest_ind(lmeanhom, l, alternative='less', equal_var=False).pvalue > 0.001:
            newgen.append('AB')
        elif np.mean(lmeanhom) > 0.5*np.mean(l):
            newgen.append('AB')
        else:
            newgen.append(row.genotype)
    if len(pd.unique(newgen)) == 1: return
    else:
        print("found: ", grp[0])
        outdf = grp[1].copy()
        outdf['newgen'] = newgen
        return outdf

with Pool(processes=64) as pool:
    newdf = pool.map(partial(getcocells, insnp=SNPPOS, fins=fins), tqdm(df3.groupby(['sample', 'bc', 'chr'])))

grp = (('MUT_11_1', 'AAACCTGAGCCCAACC', 'CUR3G'), getcelldata(df3, 'MUT_11_1', 'AAACCTGAGCCCAACC', 'CUR3G'))


if __name__ == '__main__':
    # Set constants
    INDIR = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/'
    RCDIR = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/'
    CWD = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/rtiger_out/'
    SNPFIN = {'cur':'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/strict_syn_snp_allele_readcount.depth450-650.af0.4-0.6.txt',
              'ora':'/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/ora_ref_strict_syn_snp_allele_readcount.depth450-650.af0.4-0.6.txt'}
    SYRIOUT='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out'

    SAMPLES = ('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15')
    CHRS = {'cur': ('CUR1G', 'CUR2G', 'CUR3G', 'CUR4G', 'CUR5G', 'CUR6G', 'CUR7G', 'CUR8G'),
            'ora': ('ORA1G', 'ORA2G', 'ORA3G', 'ORA4G', 'ORA5G', 'ORA6G', 'ORA7G', 'ORA8G')}

    # Read SNP markers used for CO detection
    CURSNPS = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/strict_syn_snp_allele_readcount.depth450-650.af0.4-0.6.txt', header=None)
    CURSNPS = CURSNPS.drop([3, 5, 6, 7], axis=1)
    CURSNPS.columns = ['chr', 'pos', 'ref', 'alt']
    ORASNPS = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/ora_ref_strict_syn_snp_allele_readcount.depth450-650.af0.4-0.6.txt', header=None)
    ORASNPS = ORASNPS.drop([3, 5, 6, 7], axis=1)
    ORASNPS.columns = ['chr', 'pos', 'ref', 'alt']
    SNPPOS = {grp[0]: dict(zip(grp[1]['pos'], [[0,0]]*grp[1].shape[0])) for grp in pd.concat([CURSNPS, ORASNPS]).groupby(['chr'])}

    CURLEN = {'CUR1G': 46975282,
              'CUR2G': 33806098,
              'CUR3G': 26861604,
              'CUR4G': 26096899,
              'CUR5G': 18585576,
              'CUR6G': 27136638,
              'CUR7G': 25539660,
              'CUR8G': 23045982}
    ORALEN = {'ORA1G': 47941709,
              'ORA2G': 31001264,
              'ORA3G': 27362024,
              'ORA4G': 28585890,
              'ORA5G': 18998592,
              'ORA6G': 27963823,
              'ORA7G': 24593044,
              'ORA8G': 22326149}
    GS = {'cur': sum(CURLEN.values()), 'ora': sum(ORALEN.values())}

    STEP = 100000
    WINSIZE = 1000000

    # Read the CO regions predicted by RTIGER
    curdirs = [CWD+d+"/rtiger_co_q40_lowvar/" for d in CHRS['cur']]
    with Pool(processes=4) as pool:
        curdf = pool.map(readco, curdirs)
    curdf = pd.concat(curdf)
    curdf['ref'] = 'cur'

    oradirs = [CWD+d+"/rtiger_co_q40_lowvar/" for d in CHRS['ora']]
    with Pool(processes=4) as pool:
        oradf = pool.map(readco, oradirs)
    oradf = pd.concat(oradf)
    oradf['ref'] = 'ora'

    df = pd.concat([curdf, oradf])
    df.sort_values(['sample', 'bc', 'chr'], inplace=True)
    # df.to_pickle(CWD+'df.pkl')
    df = pd.read_pickle(CWD+'df.pkl')

    # Get the coordinate of the terminal (101st SNP) marker. COs identified before this marker would be filtered out.
    # with Pool(processes=3) as pool:
    #     term = pool.map(readterm, CHRS)
    # term = pd.concat(term)
    # term.columns = ['sample', 'bc', 'chr', 'cutoff']

    # Remove chromosomes with only Het regions
    df2 = deque()
    for grp in df.groupby(['sample', 'bc', 'chr', 'ref']):
        if grp[1].shape[0] == 1:
            if grp[1]['genotype'].to_list()[0] == 'AB':
                continue
        df2.append(grp[1])
    df2 = pd.concat(df2)

    # Set all small CO regions as Het
    df2['l'] = df2['end'] - df2['start'] + 1
    df2.loc[df2['l'] < 500000, 'genotype'] = 'AB'
    df2 = mergegenotype(df2)

    # Again filter chromosomes with only Het regions
    df3 = deque()
    for grp in df2.groupby(['sample', 'bc', 'chr', 'ref']):
        if grp[1].shape[0] == 1:
            if grp[1]['genotype'].to_list()[0] == 'AB':
                continue
        df3.append(grp[1])
    df3 = pd.concat(df3)
    df3['l'] = df3['end'] - df3['start'] + 1
    # df3.to_pickle(CWD+'df3.pkl')
    df3 = pd.read_pickle(CWD+'df3.pkl')

    csnps, csnpcnt, cbinpos, cbinlen, cwindow = getwindows(CURLEN, SNPFIN['cur'], STEP=STEP, WINSIZE=WINSIZE)
    osnps, osnpcnt, obinpos, obinlen, owindow = getwindows(ORALEN, SNPFIN['ora'], STEP=STEP, WINSIZE=WINSIZE)

    df3['highhom'] = 0
    for grp in df3.groupby(['sample', 'bc', 'chr']):
        if sum(grp[1].loc[grp[1]['genotype']=='AB','l']) < 0.1 * sum(grp[1]['l']):
            df3.loc[(df3['sample']==grp[0][0]) & (df3['bc']==grp[0][1]) & (df3['chr']==grp[0][2]), 'highhom'] = 1


    fins = CWD + '{chr}/input_q40_lowvar/{chr}_{s}_{bc}_b30_q40.depth450-650.af0.4-0.6.bt2.txt'
    with Pool(processes=64) as pool:
        newdf = pool.map(partial(getcocells, insnp=SNPPOS, fins=fins), tqdm(df3.groupby(['sample', 'bc', 'chr'])))

    newdf = pd.concat(newdf)





    newdf.loc[(newdf.end - newdf.start) < 500000, 'newgen'] = 'AB'
    newdf['genotype'] = newdf['newgen']
    newdf2 = mergegenotype(newdf)
    newdf2.to_pickle(CWD+'newdf2.pkl')
    COLOR = {'AA': 'red', 'AB': 'purple', 'BB': 'blue'}
    with PdfPages(CWD + 'CO_selected_coverage_overlap_with_SNP.pdf') as pdf:
        count = 0
        fig = plt.figure(figsize=[8, 15])
        for grp in newdf2.groupby(['sample', 'bc', 'chr']):
            ax = plt.subplot2grid((30, 1), (count, 0), rowspan=2)
            acnt = {}
            with open('{INDIR}/rtiger_out/{chr}/input_q40_filter/{chr}_{s}_{bc}_b30_q40.depth350-550.af0.35-0.6.bt2.txt'.format(INDIR=INDIR, s=grp[0][0], bc=grp[0][1], chr=grp[0][2]), 'r') as fin:
                for line in fin:
                    line = line.strip().split()
                    acnt[int(line[1])] = [int(line[3]), int(line[5])]
            ax.bar(acnt.keys(), [v[0] for v in acnt.values()], width=20000, color='red')
            ax.bar(acnt.keys(), [-1*v[1] for v in acnt.values()], width=20000, color='blue')
            ax.set_ylabel('#Reads')
            ax.set_xlim([0, CHRLEN['CUR1G']])
            ax.set_title(" ".join(grp[0]))
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.set_xticks([])

            count += 2
            ax = plt.subplot2grid((30, 1), (count, 0), rowspan=2)
            ax.set_xlim([0, CHRLEN['CUR1G']])
            # ax.set_ylim([0, 1])
            ax.plot(covs[grp[0]][0].keys(), covs[grp[0]][0].values(), color='darkred', lw=0.75, label='ref cov')
            ax.plot(covs[grp[0]][1].keys(), covs[grp[0]][1].values(), color='blue', lw=0.75, label='qry cov')
            ax.plot(totcov[grp[0]].keys(), totcov[grp[0]].values(), color='black', lw=0.75, label='tot cov', ls='--')
            for row in grp[1].itertuples(index=False):
                ax.add_patch(Rectangle((row.start, 0), row.end-row.start, 1, facecolor=COLOR[row.genotype], alpha=0.5))
            ax.set_ylabel('Coverage')
            ax.hlines(avecov[grp[0]], 0, CHRLEN[grp[0][2]], label='chr mean', linestyle='dashed', color='black', lw=0.75)
            ax.hlines(cellcov[(grp[0][0], grp[0][1])], 0, CHRLEN[grp[0][2]], label='cell mean', linestyle='dotted', color='black', lw=0.75)
            # if count == 1:
            ax.legend(fontsize='small', ncol=5, loc='upper right',facecolor="none", frameon=False)

            count += 3
            if count >= 30:
                plt.subplots_adjust(left=0.1,
                                    bottom=0.02,
                                    right=0.98,
                                    top=0.98,
                                    wspace=0.4,
                                    hspace=0.5)
                # break
                pdf.savefig()
                count = 0
                fig = plt.figure(figsize=[8, 15])
        plt.tight_layout()
        plt.close()
        pdf.savefig()








    def getcoverage_SNP(grp, window_size=1000000, step_size=100000):
        """
        Get read-mapping coverage using SNP markers and separately for the two haplotypes

        Need to Optimize this. Repeatitive searching using interval tree is a bit slow
        """
        print(grp[0])
        # print(str(datetime.now()))

        snpf = CWD + '{chr}/input_q40_filter/{chr}_{sample}_{bc}_b30_q40.depth350-550.af0.35-0.6.bt2.txt'.format(sample=grp[0][0], bc= grp[0][1], chr=grp[0][2])
        rin = deque()
        qin = deque()
        selectedpos = deque()
        with open(snpf, 'r') as fin:
            for line in fin:
                line = line.strip().split()
                p = int(line[1])
                rin.append(Interval(p, p+1, int(line[3])))
                qin.append(Interval(p, p+1, int(line[5])))
                selectedpos.append(p)

        selectedpos = set(selectedpos)
        snps = SNPS.loc[SNPS['chr']==grp[0][2], 'pos']
        for snp in snps:
            if snp not in selectedpos:
                rin.append(Interval(snp, snp+1, 0))
                qin.append(Interval(snp, snp+1, 0))
        rtree = ivtree(rin)
        qtree = ivtree(qin)

        # print(str(datetime.now()))
        rmean_cov = {}
        qmean_cov = {}
        for i in range(1, CHRLEN[grp[0][2]], step_size):
            pos = rtree[i: (i+window_size)]
            try:
                rmean_cov[i+(window_size/2)] = sum([i[2] for i in pos]) / len(pos)
            except ZeroDivisionError:
                rmean_cov[i+(window_size/2)] = 0
                qmean_cov[i + (window_size / 2)] = 0
                continue
            pos = qtree[i: (i+window_size)]
            qmean_cov[i+(window_size/2)] = sum([i[2] for i in pos])/len(pos)
        # print(str(datetime.now()))
        return {grp[0]: [rmean_cov, qmean_cov]}

    with Pool(processes=80) as pool:
        covs = pool.map(getcoverage_SNP, df.groupby(['sample', 'bc', 'chr']))
    covs = {list(i.keys())[0]: list(i.values())[0] for i in covs}
    # with open(CWD+'covs.pickle.object', 'wb') as f:
    #     pickle.dump(covs, f)
    # # This covs object is saved here (/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/rtiger_out/covs.pickle.object) and can be pickle.load
    covs = pickle.load(open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/rtiger_out/covs.pickle.object', 'rb'))

    # Get Total coverage across chrs
    totcov = {}
    for k, v in covs.items():
        totcov[k] = {p: (v[0][p] + v[1][p]) for p in v[0].keys()}

    # Get average coverage for chrs
    avecov = {}
    for k, v in totcov.items():
        avecov[k] = sum(v.values())/len(v)

    # Get average coverage for cells
    cellcov = {}
    for grp in df.groupby(['sample', 'bc']):
        cov, s = 0, 0
        for c in CHRS:
            try:
                cov += avecov[(grp[0][0], grp[0][1], c)]*CHRLEN[c]
                s += CHRLEN[c]
            except KeyError:
                pass
        cellcov[(grp[0][0], grp[0][1])] = cov/s


    # Save plots for original (semi-filtered) chromosomes with homozygous regions
    COLOR = {'AA': 'red', 'AB': 'purple', 'BB': 'blue'}
    count = 1
    with PdfPages(CWD + 'CO_coverage_overlap_with_SNP.pdf') as pdf:
        fig = plt.figure(figsize=[8, 15])
        # cnt = 0
        for grp in df3.groupby(['sample', 'bc', 'chr']):
            # cnt += 1
            # if cnt ==21 : break
            # break
            ax = fig.add_subplot(10, 1, count)
            ax.set_xlim([0, CHRLEN['CUR1G']])
            # ax.set_ylim([0, 1])
            ax.plot(covs[grp[0]][0].keys(), covs[grp[0]][0].values(), color='darkred', lw=0.75, label='ref cov')
            ax.plot(covs[grp[0]][1].keys(), covs[grp[0]][1].values(), color='blue', lw=0.75, label='qry cov')
            ax.plot(totcov[grp[0]].keys(), totcov[grp[0]].values(), color='black', lw=0.75, label='tot cov', ls='--')
            for row in grp[1].itertuples(index=False):
                ax.add_patch(Rectangle((row.start, 0), row.end-row.start, 1, facecolor=COLOR[row.genotype], alpha=0.5))
            ax.set_title(" ".join(grp[0]))
            ax.set_ylabel('Coverage')
            ax.hlines(avecov[grp[0]], 0, CHRLEN[grp[0][2]], label='chr mean', linestyle='dashed', color='black', lw=0.75)
            ax.hlines(cellcov[(grp[0][0], grp[0][1])], 0, CHRLEN[grp[0][2]], label='cell mean', linestyle='dotted', color='black', lw=0.75)
            # if count == 1:
            ax.legend(fontsize='small', ncol=5, loc='upper right',facecolor="none", frameon=False)
            count += 1
            if count == 11:
                plt.tight_layout()
                pdf.savefig()
                count = 1
                fig = plt.figure(figsize=[8, 15])
        plt.tight_layout()
        pdf.savefig()
        plt.close()



    gen_cov = defaultdict(deque)
    for grp in df3.groupby(['sample', 'bc', 'chr']):
        # break
        for row in grp[1].itertuples(index=False):
            mc = 0
            cnt = 0
            for i in range(row.start, row.end, 100000):
                mc += totcov[grp[0]][((i//100000)*100000)+500001]
                cnt += 1
            gen_cov[row.genotype].append(mc/cnt)

    from scipy.stats import gaussian_kde
    fig = plt.figure()
    ax = fig.add_subplot(2, 1, 1)
    x_vals = np.linspace(0, 0.6, 200)       # Specifying the limits of our data
    density = gaussian_kde(gen_cov['AA'])
    density.covariance_factor = lambda : .5 #Smoothing parameter
    density._compute_covariance()
    ax.plot(x_vals, density(x_vals), color='red', label='AA')
    density = gaussian_kde(gen_cov['AB'])
    density.covariance_factor = lambda : .5 #Smoothing parameter
    density._compute_covariance()
    ax.plot(x_vals, density(x_vals), color='purple', label='AB')
    density = gaussian_kde(gen_cov['BB'])
    density.covariance_factor = lambda : .5 #Smoothing parameter
    density._compute_covariance()
    ax.plot(x_vals, density(x_vals), color='blue', label='BB')
    ax.legend()
    ax.set_xlabel('Coverage')
    ax.set_ylabel('Frequency')


    # Get deviation from the mean coverage of the chr
    gen_cov_dev = defaultdict(deque)
    for grp in df3.groupby(['sample', 'bc', 'chr']):
        for row in grp[1].itertuples(index=False):
            mc = 0
            cnt = 0
            for i in range(row.start, row.end, 100000):
                mc += totcov[grp[0]][((i//100000)*100000)+500001]
                cnt += 1
            gen_cov_dev[row.genotype].append((mc/cnt) - avecov[(grp[0][0], grp[0][1], grp[0][2])])

    # fig = plt.figure()
    ax = fig.add_subplot(2, 1, 2)
    x_vals = np.linspace(-0.4, 0.4, 200)       # Specifying the limits of our data
    density = gaussian_kde(gen_cov_dev['AA'])
    density.covariance_factor = lambda : .5 #Smoothing parameter
    density._compute_covariance()
    ax.plot(x_vals, density(x_vals), color='red', label='AA')
    density = gaussian_kde(gen_cov_dev['AB'])
    density.covariance_factor = lambda : .5 #Smoothing parameter
    density._compute_covariance()
    ax.plot(x_vals, density(x_vals), color='purple', label='AB')
    density = gaussian_kde(gen_cov_dev['BB'])
    density.covariance_factor = lambda : .5 #Smoothing parameter
    density._compute_covariance()
    ax.plot(x_vals, density(x_vals), color='blue', label='BB')
    ax.legend()
    ax.set_xlabel('Deviation from cell mean-coverage')
    ax.set_ylabel('Frequency')
    plt.tight_layout()
    plt.savefig(CWD+'coverage_in_het_hom_regions.pdf')

    # There are differences in the coverage of the Heterozygous and Homozygous. To remove homozygous regions that are selected only because of chromosome region dropout, I do a permutation test of the mean.
    # For a cell, I get read counts at the heterezygous positions, and then select N random positions M times to get a distribution of mean. Similarly, I get a distribution of mean for each homozygous position. If the mean read coverage of homozygous position has a different distribution than the heterozygous positions, then that positions is filtered out.

#################### TEMP FIN ############################3

################################################################################
######################### Plot CO Counts in BCs ################################
################################################################################
cocnt = deque()
for grp in df.groupby(['sample', 'bc', 'chr']):
    # cocnt.append(list(grp[0]) + [grp[1].shape[0] - 1])
    if grp[1].shape[0] == 1:
        cocnt.append(list(grp[0]) + [0])
    else:
        garb = grp[1].copy()
        garb['tmp'] = (garb['end'][:-1].reset_index(drop=True) + garb['start'][1:].reset_index(drop=True))/2
        tmpdf = garb.merge(term, on=['sample', 'bc', 'chr']).copy()
        cocnt.append(list(grp[0]) + [len(tmpdf['tmp']) - 1])
cocnt = pd.DataFrame(cocnt)
cocnt.columns = ['sample', 'bc', 'chr', 'co_count']
fig = plt.figure(figsize=[8, 10])
plti = 1
M = max(cocnt.co_count)
for grp in cocnt.groupby(['sample']):
    ax = fig.add_subplot(4, 1, plti)
    ax.set_ylim([0, 525])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    x_start = [0] * len(CHRS)
    for i in range(0, (M+1)):
        tmp = grp[1].loc[grp[1].co_count == i].copy()
        cnt = deque()
        for chr in CHRS:
            # cnt.append(garb.loc[garb.chr == chr].shape[0])
            cnt.append(tmp.loc[tmp.chr == chr].shape[0])
        ax.bar(CHRS, cnt, label=str(i), bottom=x_start, alpha=0.75)
        x_start = [x_start[j] + cnt[j] for j in range(len(cnt))]
    ax.set_xlabel(grp[0])
    ax.set_ylabel("Number of BCs")
    # ax.legend()
    plti += 1
plt.legend(loc=[0.1, -0.85], ncol=7, title="Number of COs")
plt.tight_layout()
plt.savefig(CWD + "co_counts_in_chromosomes_q40_filter.pdf")


################################################################################
################### Plot CO Frequency along chromosome #########################
################################################################################


coloc = pd.DataFrame()
for grp in df.groupby(['sample', 'bc', 'chr']):
    if grp[1].shape[0] == 1: continue
    loc = (grp[1]['end'][:-1].reset_index(drop=True) + grp[1]['start'][1:].reset_index(drop=True))/2
    loc = pd.DataFrame(loc)
    loc.columns = ['pos']
    loc['sample'] = grp[0][0]
    loc['bc'] = grp[0][1]
    loc['chr'] = grp[0][2]
    loc = loc[['sample', 'bc', 'chr', 'pos']]
    coloc = pd.concat([coloc, loc])

fig = plt.figure(figsize=[12, 9])
plti = 1
linestyle = {'WT_1': '-k', 'WT_19': '-.k', 'MUT_11_1': '-r', 'MUT_15': '-.r'}
for grp in coloc.groupby('chr'):
    ax = fig.add_subplot(4, 2, plti)
    ax.set_ylim([0, 25])
    for grp2 in grp[1].groupby('sample'):
        ax.set_xlim([0, CHRLEN['CUR1G']])
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        bins = list(range(1, CHRLEN[grp[0]], 1000000)) + [CHRLEN[grp[0]]]
        nbin = len(bins) - 1
        poscnt = pd.cut(sorted(grp2[1]['pos']), bins).value_counts().values
        x = (np.array(bins)[0:-1] + np.array(bins)[1:])/2
        ax.plot(x, poscnt, linestyle[grp2[0]], label=grp2[0])
        ax.set_ylabel('#CO')
    ax.set_xlabel(grp[0])
    ax.legend(fontsize='x-small', loc=1)
    plti += 1
plt.tight_layout()
plt.savefig(CWD +'co_counts_along_chromosome_q40_filter.pdf')



################################################################################
############# Plot Chromosome CO-Count vs Average read depth ###################
################################################################################


read_depth = pd.DataFrame()
for sample in SAMPLES:
    sdf = pd.DataFrame()
    bcs = os.listdir("{}/{}/barcodes/".format(INDIR, sample))
    for bc in bcs:
        bcdf = deque()
        with open("{}/{}/barcodes/{}/mean_read_cov_q30_Q40.txt".format(INDIR, sample, bc), 'r') as fin:
            for line in fin:
                if 'CUR' in line:
                    line = line.strip().split()
                    bcdf.append([sample, bc, line[0], round(float(line[1]), 2)])
        bcdf = pd.DataFrame(bcdf)
        sdf = pd.concat([sdf, bcdf])
    read_depth = pd.concat([read_depth, sdf])
read_depth.columns = ['sample', 'bc', 'chr', 'depth']

cord = cocnt.merge(read_depth, on=['sample', 'bc', 'chr'])
import seaborn as sns
sns.catplot(x='co_count', y='depth', data=cord, hue='sample',
            col='chr', col_wrap=2, kind='violin', inner=None,
            height=3, aspect=2.5, linewidth=0)
plt.tight_layout()
plt.tight_layout()
plt.savefig(CWD+'CO_read_depth_distribution_q40_filter.pdf')


################################################################################
######################## Find co-occuring COs ##################################
################################################################################


from matplotlib.ticker import AutoMinorLocator
hiloc = deque()
for grp in coloc.groupby('chr'):
    tmp = grp[1].copy()
    pos = tmp['pos'].to_list()
    for p in pos:
        s = sum(abs(tmp['pos'] - p) < 20000)
        if s >= 1:
            hiloc.append([grp[0], p, s])
hiloc = pd.DataFrame(hiloc)
hiloc.columns = ['chr', 'pos', 'co']
fig = plt.figure(figsize=[25, 20])
plti = 1
for grp in hiloc.groupby('chr'):
    ax = fig.add_subplot(8, 1, plti)
    ax.set_ylim([0, 25])
    ax.set_xlim([0, CHRLEN['CUR1G']])
    ax.set_title(grp[0])
    ax.set_axisbelow(True)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(which='both', axis='both')
    ax.vlines(x=grp[1]['pos'], ymin=0, ymax=grp[1]['co'], linewidth=0.5, color='grey', zorder=1)
    ax.scatter(grp[1]['pos'], grp[1]['co'], color='black', s=4, zorder=3)
    ax.set_ylabel("CO Count")
    snps = SNPS.loc[SNPS['chr'] == grp[0]]
    ax.vlines(x=snps['pos'], ymin=0, ymax=25, color='lightblue', linewidth=0.1, zorder=0, alpha=0.1)
    plti += 1
ax.set_xlabel("Chromosome position (in Mbp)")
plt.tight_layout()
plt.savefig(CWD+'CO_clusters_q40_filter.pdf')
plt.savefig(CWD+'CO_clusters_q40_filter.png')
plt.close()

## Select plots for CO Clusters
hiloc2 = hiloc.loc[hiloc['co'] > 5].copy()
hiloc2.sort_values(['chr', 'pos'], inplace=True)
hiloc2.reset_index(inplace=True, drop=True)
posloc = [41]

c = hiloc2.iat[41, 0]
p = hiloc2.iat[41, 1]
fs = coloc2.loc[(coloc2['pos'] >= p-20000) & (coloc2['pos'] < p+20000) & (coloc2['chr'] == c)]
fins = ['{CWD}/{c}/rtiger_co/{c}_{s}_{bc}_b30_q10.bt2.txt/GenotypePlot_{c}_{s}_{bc}_b30_q10.bt2.txt.pdf'.format(CWD=CWD, c=c, s=row[0], bc=row[1]) for row in fs.itertuples(index=False)]
fout = '/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/tests/temp.pdf'
mergepdf(fins, fout)


#############################################################################
############### Chromosomes with large Homozygous regions ###################
#############################################################################


def gethighhom(df, CWD, hetrat=0.5):
    fins = deque()
    cnt = 0
    for grp in df.groupby(['sample', 'bc', 'chr']):
        cnt += 1
        tl = sum(grp[1]['end'] - grp[1]['start'] + 1)
        het = grp[1].loc[grp[1]['genotype']=='AB']
        hetl = sum(het['end'] - het['start'] + 1)
        if hetl/tl < hetrat:
            fins.append('{CWD}/{c}/rtiger_co_q40_filter/{c}_{s}_{bc}_b30_q40.depth350-550.af0.35-0.6.bt2.txt/GenotypePlot_{c}_{s}_{bc}_b30_q40.depth350-550.af0.35-0.6.bt2.txt.pdf'.format(CWD=CWD, c=grp[0][2], s=grp[0][0], bc=grp[0][1]))
    return (fins, cnt)

fins, cnt = gethighhom(df, CWD, 0.5)
# Number of Chromosomes tested: 11378 of which 526 (4.6%) have > 50% homozygous regions
mergepdf(fins, CWD+'high_homo_chrs_q40_hom50.pdf')

fins, cnt = gethighhom(df, CWD, 0.1)
# Number of Chromosomes tested: 11378 of which 168 (1.4%) have > 90% homozygous regions
mergepdf(fins, CWD+'high_homo_chrs_q40_filter_hom90.pdf')
alls = [grp[0][0] for grp in df.groupby(['sample','bc', 'chr'])]
selecteds = [s for f in fins for s in SAMPLES if s+"_" in f]

allchrs = [grp[0][2] for grp in df.groupby(['sample','bc', 'chr'])]
selectedchrs = [c for f in fins for c in CHRS if c in f]

fig = plt.figure(figsize=[12,8])
ax = fig.add_subplot(2, 2, 1)
ax.bar(SAMPLES, [Counter(alls)[k] for k in SAMPLES], label='Total')
ax.bar(SAMPLES, [Counter(selecteds)[k] for k in SAMPLES], label='Homo')
ax.set_ylabel('Chr Count')
ax.legend()
ax = fig.add_subplot(2, 2, 2)
ax.bar(SAMPLES, [Counter(selecteds)[k]*100/Counter(alls)[k] for k in SAMPLES], label='Total')
ax.set_ylabel('Percent Homo')

ax = fig.add_subplot(2, 2, 3)
ax.bar(CHRS, [Counter(allchrs)[k] for k in CHRS], label='Total')
ax.bar(CHRS, [Counter(selectedchrs)[k] for k in CHRS], label='Homo')
ax.set_ylabel('Chr Count')
ax.legend()
ax = fig.add_subplot(2, 2, 4)
ax.bar(CHRS, [Counter(selectedchrs)[k]*100/Counter(allchrs)[k] for k in CHRS], label='Total')
ax.set_ylabel('Percent Homo')
plt.tight_layout()
plt.savefig(CWD+'homo_dist_sample_chr.pdf')


################################################################################
############################## Plot CO Density #################################
################################################################################


fig = plt.figure(figsize=[8, 16])
plti = 1
for grp in coloc.groupby('chr'):
    ax = fig.add_subplot(8, 1, plti)
    data = sorted(grp[1]['pos'])
    bins = list(range(1, CHRLEN[grp[0]], 1000000)) + [CHRLEN[grp[0]]]
    y = pd.cut(data, bins).value_counts().values
    x = [(bins[i] + bins[i+1])/2 for i in range(len(bins) - 1)]
    ax.plot(x, y)

    snps = SNPS.loc[SNPS['chr'] == grp[0]]
    ax.vlines(x=snps['pos'], ymin=0, ymax=max(y), color='lightblue', linewidth=0.1,zorder=0, alpha=0.1)
    ax.set_xlim([0, CHRLEN['CUR1G']])
    ax.set_title(grp[0], size=10)
    ax.set_ylabel('#CO')
    plti+=1


# # Read CUR6G reversed COs
# os.chdir(CWD+'/CUR6G/reversed_rtiger_co_q40/')
# dirs = [d for d in os.listdir() if os.path.isdir(d)]
# revdf = pd.DataFrame()
# for d in dirs:
#     s = [s for s in SAMPLES if s+'_' in d][0]
#     # print(s)
#     bc = re.findall(r'[a-zA-Z]{16}', d)[0]
#     g = [f for f in os.listdir(d) if 'CompleteBlock-state' in f][0]
#     bcdf = pd.read_table(d + '/' + g, header=None)
#     bcdf.columns = ['chr', 'start', 'end', 'genotype']
#     bcdf['sample'] = s
#     bcdf['bc'] = bc
#     bcdf = bcdf[['sample', 'bc', 'chr', 'start', 'end', 'genotype']]
#     revdf = pd.concat([revdf, bcdf])
#
# revcoloc = pd.DataFrame()
# for grp in revdf.groupby(['sample', 'bc', 'chr']):
#     if grp[1].shape[0] == 1: continue
#     loc = (grp[1]['end'][:-1].reset_index(drop=True) + grp[1]['start'][1:].reset_index(drop=True))/2
#     loc = pd.DataFrame(loc)
#     loc.columns = ['pos']
#     loc['sample'] = grp[0][0]
#     loc['bc'] = grp[0][1]
#     loc['chr'] = grp[0][2]
#     loc = loc[['sample', 'bc', 'chr', 'pos']]
#     revcoloc = pd.concat([revcoloc, loc])
#
#
# ax = fig.add_subplot(9, 1, plti)
# data = sorted(revcoloc['pos'])
# bins = list(range(1, CHRLEN['CUR6G'], 1000000)) + [CHRLEN['CUR6G']]
# y = pd.cut(data, bins).value_counts().values
# x = [(bins[i] + bins[i+1])/2 for i in range(len(bins) - 1)]
# ax.plot(x, y)
# ax.set_xlim([0, CHRLEN['CUR1G']])
# ax.set_title('Reversed CUR6G', size=10)
ax.set_ylabel('#CO')
ax.set_xlabel('Chromosome Postion')
plt.tight_layout()

plt.savefig(CWD + '/CO_density_filter_q40.pdf')


############################################################################
############### Filter out noisy low-confidence COs ########################
############################################################################


################################################################################
########## Filter homozygous regions based on noise from other parent ##########
################################################################################



def checkhomo(indf):
    outdf = pd.DataFrame()
    refratio = deque()
    qryratio = deque()
    homrat = deque()
    for grp in indf.groupby(['sample', 'bc', 'chr']):
        # break
        tmp = deque()
        with open(CWD + '/{chr}/input_q40_filter/{chr}_{s}_{bc}_b30_q40.depth350-550.af0.35-0.6.bt2.txt'.format(chr=grp[0][2], s=grp[0][0], bc=grp[0][1])) as fin:
           for row in grp[1].itertuples(index=False):
               refcnt = 0
               qrycnt = 0
               # if row.genotype == 'AB':
               #     homrat = tmp.append(list(row))
               for line in fin:
                   line = line.strip().split()
                   p = int(line[1])
                   if p < row.start: continue
                   if p > row.end: break
                   refcnt += int(line[3])
                   qrycnt += int(line[5])
                homrat.append(refcnt/(refcnt+qrycnt))
               # if row.genotype == 'AA':
               #     refratio.append(refcnt/(refcnt+qrycnt))
               #     if refcnt/(refcnt+qrycnt) > 0.99: tmp.append(list(row))
               #     else:
               #         t = list(row)
               #         t[5] = 'AB'
               #         tmp.append(t)
               # elif row.genotype == 'BB':
               #     qryratio.append(qrycnt/(refcnt+qrycnt))
               #     if qrycnt/(refcnt+qrycnt) > 0.99: tmp.append(list(row))
               #     else:
               #         t = list(row)
               #         t[5] = 'AB'
               #         tmp.append(t)
        # tmp = pd.DataFrame(tmp)
        # tmp.columns = row._fields
        # outdf = pd.concat([outdf, tmp])

# Plot was generated before changing checkhomo to calculate ratio for all positions
ratvau, bas = np.histogram(list(refratio) + list(qryratio), bins=[r/100 for r in range(70, 101)])
ratcum = sum(ratvau) - np.cumsum(ratvau)
plt.bar(bas[:-1]+0.005, ratvau, width=0.01)
plt.xlabel('Homozygous allele ratio')
plt.ylabel('Number of homozygous regions')
plt.grid(which='both', axis='x')
plt.twinx()
plt.plot(bas[1:], ratcum, color='black', label='cumulative')
plt.ylabel('Cumulative number of homozygous regions')
plt.legend()
plt.grid(which='both', axis='y')
plt.savefig(CWD+'good_cov_CO_hom_ratio.pdf')

#
indf['homrat'] = homrat
bins = [r/50 for r in range(51)]
plt.hist(indf.loc[indf['genotype']=='AB', 'homrat'], bins, alpha=0.5, color='black')
plt.hist(indf.loc[indf['genotype']=='AA', 'homrat'], bins, alpha=0.5, color='red')
plt.hist(indf.loc[indf['genotype']=='BB', 'homrat'], bins, alpha=0.5, color='blue')







fins = deque()
for row in outdf.itertuples(index=False):
    fins.append(CWD + '{chr}/rtiger_co_q40_filter/{chr}_{sample}_{bc}_b30_q40.depth350-550.af0.35-0.6.bt2.txt/GenotypePlot_{chr}_{sample}_{bc}_b30_q40.depth350-550.af0.35-0.6.bt2.txt.pdf'.format(sample=row[0], chr=row[2], bc=row[1]))
    break
CUR6G/rtiger_co_q40_filter/CUR6G_MUT_15_GTCTTCGAGTGTCCAT_b30_q40.depth350-550.af0.35-0.6.bt2.txt/GenotypePlot_CUR6G_MUT_15_GTCTTCGAGTGTCCAT_b30_q40.depth350-550.af0.35-0.6.bt2.txt.pdf

plt.hist(refratio, bins=[i/50 for i in range(0, 50)], label='Parent 1', color='blue', alpha=0.5, histtype='bar')
plt.hist(qryratio, bins=[i/50 for i in range(0, 50)], label='Parent 2', color='red', alpha=0.5, histtype='bar')
plt.xlabel('Ratio of reads corresponding to parent of the homozygous region')
plt.ylabel('Number of Regions')
plt.legend()
plt.savefig(CWD+'read_ratio_in_homozygous_reagions_filter.pdf')
plt.savefig(CWD+'read_ratio_in_homozygous_reagions_filter.png')

with open(CWD+'/homo_regions_ratio.txt', 'w') as fout:
    fout.write('Parent_1_ratio: ')
    for r in refratio: fout.write(str(round(r, 6)) + ',')
    fout.write('\nParent_2_ratio: ')
    for r in qryratio: fout.write(str(round(r, 6)) + ',')



# If small Het region is surrounded by same Hom region, then set Het -> Hom
df3 = pd.DataFrame()
for grp in df2.groupby(['sample', 'bc', 'chr']):
    tmp = grp[1].copy()
    if tmp.shape[0] == 1:
        df3 = pd.concat([df3, tmp])
        continue
    g = tmp['genotype'].to_numpy()
    heti = np.where(g=='AB')[0]
    s = tmp['l'].to_numpy()
    n = len(g)-1
    try:
        for i in heti:
            if s[i] > 500000: continue
            if i == 0: g[i] = g[i+1]
            elif i == n: g[i] = g[i-1]
            elif g[i-1] == g[i+1]:
                if g[i-1] == 'AB':
                    print('ERROR')
                g[i] = g[i-1]
    except:
        print(grp[1])
        break
    tmp['genotype'] = g
    df3 = pd.concat([df3, tmp])
df3 = mergegenotype(df3)
df3['l'] = df3['end']- df3['start'] + 1

garb = mergegenotype(df3)
garb['l'] = garb['end']- garb['start'] + 1




for row in df3.itertuples(index=False):
    break
sdf = pd.read_table('{CWD}{chr}/input_q40/{chr}_{sample}_{bc}_b30_q40.bt2.txt'.format(CWD=CWD, chr=row.chr, sample=row.sample, bc=row.bc), header=None)
sdf.columns = ['chr', 'pos', 'ref', 'ra', 'qry', 'qa']
garb = sdf.loc[(sdf['pos']>row.start) & (sdf['pos']<row.end)]




######## Plot CO Size Distribution #################

fig = plt.figure(figsize=[10, 8])
plti = 1
for grp in df.groupby(['sample']):
    ax = fig.add_subplot(4, 1, plti)
    for chr in CHRS:
        bad_chr = cocnt.loc[(cocnt['sample']==grp[0]) & (cocnt['chr']==chr) & (cocnt['co_count']==0)]
        tmp = grp[1].loc[grp[1]['chr'] == chr]
        tmp = tmp.loc[~(tmp['bc'].isin(bad_chr['bc']))]
        l = tmp['end'] - tmp['start'] + 1
        y, x, _ = ax.hist(l, bins=range(1, CHRLEN[chr], 500000), histtype='stepfilled', alpha=0, label=chr)
        ax.plot([j + 250000 for j in x[:-1]], y)
    ax.set_title(grp[0])
    ax.set_xlabel('Genotype Size')
    ax.set_ylabel("Frequency")

    # break
    plti += 1
plt.legend(loc=[0.1, -0.85], ncol=4, title="Chromosome")
plt.tight_layout()



################################################################################
############## Read coverage along chromosome for all samples ##################
################################################################################
snpsd = defaultdict(dict)
for row in SNPS.itertuples(index=False):
    snpsd[row[0]][row[1]] = 0

d = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/'
for sample in SAMPLES:
    ssnp = snpsd.copy()
    dirs = sorted([dir for dir in os.listdir(d + sample) if os.path.isdir(d+sample+'/'+dir)])
    for bc in dirs:
        print(bc)
        with open(d+sample+'/'+bc+'/{sample}_{bc}_b30_q40.bt2.txt'.format(sample=sample, bc=bc), 'r') as f:
            for line in f:
                line = line.strip().split()
                ssnp[line[0]][int(line[1])] += int(line[3]) + int(line[5])

    fig = plt.figure()
    plti = 1
    for chr in CHRS:
        ax = fig.add_subplot(8, 1, plti)
        pos = np.array(list(ssnp[chr].keys()), dtype='int')
        bins = list(range(1, CHRLEN[chr], 10000)) + [CHRLEN[chr]]
        bcnt = {}
        for i in range(len(bins) -1):
            b = bins[i]
            xs = np.argwhere((pos>b) & (pos < (b+100000)))
            s = 0
            for x in xs:
                s += ssnp[chr][int(pos[x])]
            try:
                bcnt[b] = s/(len(xs))
            except ZeroDivisionError as e:
                bcnt[b] = 0

        # xbin = pd.cut(list(ssnp[chr].keys()), bins)
        # x = sorted(np.random.choice(list(ssnp[chr].keys()), 10000))
        keys = sorted(list(bcnt.keys()))
        ax.plot(keys, [bcnt[k] for k in keys])
        ax.set_xlim([0, CHRLEN['CUR1G']])
        plti += 1

    plt.savefig(d+'rtiger_out/sequenced_co_marker_density_{}.pdf'.format(sample))

rsnp = defaultdict(dict)
qsnp = defaultdict(dict)
for row in SNPS.itertuples(index=False):
    rsnp[row[0]][row[1]] = 0
    qsnp[row[0]][row[1]] = 0


d = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/'
fig = plt.figure()
plti = 1
qsnp0 = deque()
for sample in SAMPLES:
    rsnp = defaultdict(dict)
    qsnp = defaultdict(dict)
    for row in SNPS.itertuples(index=False):
        rsnp[row[0]][row[1]] = 0
        qsnp[row[0]][row[1]] = 0

    dirs = sorted([dir for dir in os.listdir(d + sample) if os.path.isdir(d+sample+'/'+dir)])[1:100]
    for bc in dirs:
        print(bc)
        with open(d+sample+'/'+bc+'/{sample}_{bc}_b30_q40.depth350-550.af0.35-0.6.bt2.txt'.format(sample=sample, bc=bc), 'r') as f:
            for line in f:
                line = line.strip().split()
                rsnp[line[0]][int(line[1])] += int(line[3])
                qsnp[line[0]][int(line[1])] += int(line[5])

    for chr, poss in qsnp.items():
        for pos, v in poss.items():
            if v == 0:
                qsnp0.append(chr + str(pos))
    for chr in CHRS:
        ax = fig.add_subplot(4, 8, plti)
        max_x = max(max(rsnp[chr].values()), max(qsnp[chr].values()))
        ax.hist(rsnp[chr].values(), bins=range(max_x), alpha=0.5, label='P1')
        ax.hist(qsnp[chr].values(), bins=range(max_x), alpha=0.5, label='P2')
        if sample == 'MUT_15': ax.set_xlabel('read count')
        if chr == 'CUR1G': ax.set_ylabel('{}\n#Markers'.format(sample))
        if sample == 'WT_1': ax.set_title(chr)
        plti += 1