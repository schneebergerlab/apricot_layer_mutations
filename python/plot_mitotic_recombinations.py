#!/usr/bin/env python3
import numpy as np
# import pandas as pd
from matplotlib import pyplot as plt
import os, sys
from collections import defaultdict, Counter, deque
from pathlib import Path
from scipy.stats import binom_test
from matplotlib.backends.backend_pdf import PdfPages
from multiprocessing import Pool
from functools import partial

SAMPLES     = ("WT_1", "WT_19", "MUT_11_1", "MUT_15")
CWD         = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/'
BASEDICT    = {'A': 4, 'C': 5, 'G': 6, 'T': 7}


def ax_rem_spine(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.grid(which='major', axis='x', linestyle='--')
    ax.set_axisbelow(True)

def getsampleplot(sample, CWD, snppos, binspos, binsmid):
    count = 0
    os.chdir(CWD + sample)
    # print("Analysing: ", os.getcwd())
    # return
    with PdfPages('mitotic_recombination.pdf') as pdf:
        for bc in [d for d in Path('.').iterdir() if d.is_dir()]:
            # if str(bc)[0] != 'A':
            #     continue
            # print(bc)
            finname = [f for f in os.listdir("{}/".format(bc)) if 'read_counts_b30_q40' in f][0]
            with open(str(bc) + os.sep + finname, 'r') as fin:
                bcdata = defaultdict(dict)
                for line in fin:
                    line = line.strip().split()
                    if line[3] == '0': continue
                    try:
                        bcdata[line[0]][int(line[1])] = (int(line[BASEDICT[snppos[line[0]][int(line[1])][0]]]), int(line[BASEDICT[snppos[line[0]][int(line[1])][1]]]))
                    except KeyError:
                        pass

            bcpred = defaultdict(dict)
            for chr, indices in binspos.items():
                state = deque()
                ratio = deque()
                acnt = deque()
                bcnt = deque()
                poslen = deque()
                for index, pos in indices.items():
                    a = 0
                    b = 0
                    pl = 0
                    for p in pos:
                        try:
                            val = bcdata[chr][p]
                            a += val[0]
                            b += val[1]
                            pl += 1
                        except KeyError as e:
                            continue

                    acnt.append(a)
                    bcnt.append(b)
                    poslen.append(pl)
                    if a + b == 0:
                        state.append(0)
                        ratio.append(0)
                        continue
                    ratio.append((2*(a/(a+b))) - 1)
                    if a+b < 100:
                        state.append(0)
                        continue
                    if a/(a+b) > 0.125 and a/(a+b) < 0.875:
                        state.append(0)
                        continue
                    pbin = binom_test(a, a+b, 0.5, 'two-sided')
                    if pbin > 0.01: state.append(0)
                    elif a > b: state.append(1)
                    else: state.append(-1)

                rcs = [acnt[i] + bcnt[i] for i in range(len(acnt)) if (acnt[i] + bcnt[i]) > 0]
                if len(rcs) > 0:
                    rccut = 0.2*(sorted(rcs)[int(len(rcs)/2)])
                    for i in range(len(state)):
                        if acnt[i] + bcnt[i] < rccut:
                            state[i] = 0
                if '111' in ''.join(map(str, state)) or '-1-1-1' in ''.join(map(str, state)):
                    # count += 1
                    # if count == 100:
                    #     return
                    print('Found Mitotic Recombination: {} {} {}'.format(sample, str(bc), chr))
                    fig = plt.figure(figsize=[12, 8])
                    xpos = list(binsmid[chr].values())

                    ax = fig.add_subplot(4, 1, 1)
                    ax_rem_spine(ax)
                    ax.set_ylim([-1.01, 1.01])
                    ax.set_xlim([0, max(xpos)+10000])
                    ax.set_title(str(bc) + ' ' + chr)
                    refy = deque()
                    alty = deque()
                    for i in range(len(acnt)):
                        try:
                            refy.append(acnt[i]/poslen[i])
                            alty.append(-bcnt[i]/poslen[i])
                        except ZeroDivisionError:
                            refy.append(0)
                            alty.append(0)
                    ax.bar(xpos, refy, color='red', width=100000)
                    ax.bar(xpos, alty, color='blue', width=100000)
                    ax.set_xticklabels([str(int(i)) for i in ax.get_xticks()])
                    ax.set_ylabel('reads per sequenced\nmarker in window')

                    ax = fig.add_subplot(4, 1, 2)
                    ax_rem_spine(ax)
                    maxy = max(max(acnt), max(bcnt))
                    ax.set_ylim([-maxy, maxy])
                    ax.set_xlim([0, max(xpos)+10000])
                    refy = deque()
                    alty = deque()
                    for i in range(len(acnt)):
                        refy.append(acnt[i])
                        alty.append(-bcnt[i])
                    ax.bar(xpos, refy, color='red', width=100000)
                    ax.bar(xpos, alty, color='blue', width=100000)
                    ax.set_xticklabels([str(int(i)) for i in ax.get_xticks()])
                    ax.set_ylabel('total reads \n in window ')

                    ax = fig.add_subplot(4, 1, 3)
                    ax_rem_spine(ax)
                    ax.set_ylim([-1.01, 1.01])
                    ax.set_xlim([0, max(xpos)+10000])
                    refy = np.ma.masked_where(np.array(ratio) < 0.75, np.array(ratio))
                    hety = np.ma.masked_where((np.array(ratio) < -0.75) | (np.array(ratio) > 0.75) , np.array(ratio))
                    alty = np.ma.masked_where(np.array(ratio) > -0.75, np.array(ratio))
                    ax.hlines([-0.75, 0.75], 0, max(xpos)+10000, linestyles='dashed', colors='grey')
                    ax.plot(xpos, refy, 'red', xpos, alty, 'blue', xpos, hety, 'purple')
                    ax.set_xticklabels([str(int(i)) for i in ax.get_xticks()])
                    ax.set_ylabel('Ratio of CUR \n to ORA reads')

                    ax = fig.add_subplot(4, 1, 4)
                    ax_rem_spine(ax)
                    ax.set_ylim([-1.01, 1.01])
                    ax.set_xlim([0, max(xpos)+10000])
                    ax.plot(xpos, state)
                    ax.set_xticklabels([str(int(i)) for i in ax.get_xticks()])
                    ax.set_ylabel('Genomic state')

                    plt.tight_layout()
                    plt.subplots_adjust(hspace=0.5)
                    pdf.savefig()
                    plt.close()
                    # return


snppos = defaultdict(dict)
# with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/syn_snp.txt', 'r') as fin:
with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/strict_syn_snp.selected.txt', 'r') as fin:
    for line in fin:
        line = line.strip().split()
        if line[3] == 'N' or line[4] == 'N': continue
        snppos[line[0]][int(line[1])] = (line[3], line[4])

# print('Read SNP positions')

WS = 1000000 # Window Size
SS = 100000 # Step size
binspos = defaultdict(dict)
binsmid = defaultdict(dict)
for chr, pos in snppos.items():
    maxpos = max(list(pos.keys()))
    bins = list(zip(range(1, maxpos, SS), range(WS, maxpos, SS)))
    bp = {}
    bm = {}
    for i in range(len(bins)):
        bp[i] = [p for p in pos.keys() if p>=bins[i][0] and p<=bins[i][1]]
        bm[i] = (bins[i][0] + bins[i][1])/2
    binspos[chr] = bp
    binsmid[chr] = bm

# print('Generated Bins')

os.chdir(CWD)


for sample in SAMPLES:
    getsampleplot(sample, CWD=CWD, snppos=snppos, binspos=binspos, binsmid=binsmid)


# with Pool(processes=4) as pool:
#     out = pool.map(partial(getsampleplot, CWD=CWD, snppos=snppos, binspos=binspos, binsmid=binsmid), SAMPLES)
#
#
# with PdfPages('test.pdf') as pdf:
#     for i in range(200):
#         plt.hist([1,2,3,4,5])
#         pdf.savefig()
#         plt.close()