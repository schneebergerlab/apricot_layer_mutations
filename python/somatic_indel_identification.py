import argparse

def write_pos_only_indels(fin, fout, CAN_N, RC_MIN, RC_MAX, CAN_CNT):
    '''
    Select all positions that one alternate substitution allele (SNPs).
    '''
    from collections import deque

    with open(fin, 'r') as f:
        count = 0
        df = deque()
        for line in f:
            count += 1
            if count%1000000 == 0: print(count)
            line = line.strip().split()
            if len(line) < 11: continue
            if int(line[3]) < RC_MIN: continue
            if int(line[3]) > RC_MAX: continue

            # Filter out positions where non-reference bases are also present, as these are noisy positions and we would not select indels there.
            hf = 0
            for i in [4, 5, 6, 7]:
                if int(line[i]) >= CAN_N: hf += 1
            if hf > 1: continue

            # Filter out positions that do not have exactly one indel with >= CAN_N reads
            hf = 0
            for i in range(10, len(line), 2):
                if int(line[i]) >= CAN_N: hf += 1
            if hf != 1: continue

            try:
                # Get indel type, alternate allele readcount and alternate allele ratio
                indel = ''
                rc = 0
                ar = 0
                for i in range(10, len(line), 2):
                    if int(line[i]) >= CAN_N:
                        indel = line[i-1]
                        rc = int(line[i])
                        ar = int(line[i])/int(line[3])

                if rc == 0 or ar == 0 or indel == '': continue
                df.append([line[0], line[1], line[1], line[2], indel, rc, ar])

            except KeyError as e:
                pass
            # if line[1] == '9163':
            #     print(rc, ar, indel)
            #     break

    import pandas as pd
    df = pd.DataFrame(df)
    df[1] = df[1].astype(int)
    df[2] = df[2].astype(int)
    df[5] = df[5].astype(int)
    df_cand = df.iloc[0:CAN_CNT] if CAN_CNT != -1 else df.copy()
    df_cand.to_csv(fout+'.regions', sep=' ', index=False, header=False)
    df_bed = pd.DataFrame({'chr'  : df_cand[0],
                           'start': df_cand[1]-1,
                           'end'  : df_cand[2],
                           'ref'  : df_cand[3],
                           'alt'  : df_cand[4],
                           'vaf'  : df_cand[5],
                           'var'  : df_cand[6]})
    # df_bed.to_csv(foutname+'.bed', sep='\t', index=False, header=False)
    df_bed.sort_values(['chr', 'start', 'end'], inplace=True)
    df_bed.to_csv(fout +'.sorted.bed', sep='\t', index=False, header=False)


def getcandidate(args):
    import sys
    if len(args.f) == 1:
        sys.exit('Need at least two samples to find somatic mutations')

    if len(args.s) != 1:
        if len(args.s) != len(args.f):
            sys.exit('Prefix for output file should either be one string or each input file should get one prefix')

    FINS = [i.name for i in args.f]
    CAN_N = args.n
    RC_MIN = args.m
    RC_MAX = args.M
    PRES = [args.s[0]+'_' for i in range(len(FINS))] if len(args.s) == 1 else [i+'_' for i in args.s]
    CAN_CNT = args.N
    N = len(FINS)
    CORES = args.C


    from collections import deque, defaultdict
    import pandas as pd
    from matplotlib import pyplot as plt
    from multiprocessing import Pool
    from functools import partial
    # from upsetplot import UpSet


    FINS = ['WT_1/bam_read_counts_b30_q10.bt2.txt', 'WT_19/bam_read_counts_b30_q10.bt2.txt','MUT_11_1/bam_read_counts_b30_q10.bt2.txt','MUT_15/bam_read_counts_b30_q10.bt2.txt']
    PRES=['WT_1_bt2_', 'WT_19_bt2_', 'MUT_11_1_bt2_', 'MUT_15_bt2_']

    fouts = deque()
    for i in range(N):
        # For each sample select and write only SNP positions
        if '/'.join(FINS[i].split("/")[:-1]) == '':
            fout = './' + "/{}".format(PRES[i]) + 'indel_candidate'
        else:
            fout = '/'.join(FINS[i].split("/")[:-1]) + "/{}".format(PRES[i]) + 'indel_candidate'
        fouts.append(fout)

    with Pool(processes=CORES) as pool:
        pool.starmap(partial(write_pos_only_indels, CAN_N=CAN_N, RC_MIN=RC_MIN, RC_MAX=RC_MAX, CAN_CNT=CAN_CNT), zip(FINS, fouts))


    ## Reading Bed files generated above to perform filtersing
    import pybedtools
    beds = deque()
    for f in fouts:
        beds.append(pybedtools.BedTool(f + '.sorted.bed'))


    ## Select positions that are mutated in all samples
    beds_reg = defaultdict(int)
    for bed in beds:
        for b in bed:
            beds_reg['_'.join(list(b)[:5])] += 1
    bad_pos = set([k for k, v in beds_reg.items() if v == N])


    ## Output positions that are not mutated in all samples
    for i in range(N):
        with open(fouts[i]+'.sorted.filtered.bed', 'w') as f:
            for b in beds[i]:
                if '_'.join(list(b)[:5]) not in bad_pos:
                    f.write('\t'.join(list(b)) + '\n')


def multicell(args):
    F1 = args.pos.name
    F3 = args.rc_file_list.name
    M  = args.m
    PRES  = args.o

    import pandas as pd
    from collections import Counter

    f1 = pd.read_table(F1, header=None)
    f3 = pd.read_table(F3, header=None)
    f1.columns = ['chr', 'pos0', 'pos', 'ref', 'alt', 'vac', 'var']
    f3.columns = ['fin']

    f1_alt = {row.chr + '_' + str(row.pos): [row.alt, 0] for row in f1.itertuples(index=False)}

    loci = set(f1_alt.keys())
    for bc in f3.fin:
        with open(bc) as fin:
            for line in fin:
                line = line.strip().split()
                if line[0] + '_' + line[1] not in loci: continue
                if len(line) < 10: continue
                for i in range(9, len(line), 2):
                    if line[i] == f1_alt[line[0] + '_' + line[1]][0]:
                        f1_alt[line[0] + '_' + line[1]][1] += 1

    f1_loci = list(f1['chr'].astype(str) + '_' + f1['pos'].astype(str))
    f1_g = [True if f1_alt[i][1] >= M else False for i in f1_loci]
    f1_g_can = f1.loc[f1_g].copy()
    if '/'.join(F1.split("/")[:-1]) == '':
        fout = './' + "/{}".format(PRES) + '_' + F1.split("/")[-1]
    else:
        fout = '/'.join(F1.split("/")[:-1]) + "/{}".format(PRES) + '_' + F1.split("/")[-1]
    f1_g_can.to_csv(path_or_buf=fout,
                    sep='\t',
                    header=False,
                    index=False)

    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=[5, 3])
    ax = fig.add_subplot(1, 1, 1)
    ax.bar(x= list(Counter([i[1] for i in f1_alt.values()]).keys()),
           height=list(Counter([i[1] for i in f1_alt.values()]).values()))
    ax.set_yscale('log')
    ax.set_xlabel('Number of cells')
    ax.set_ylabel('Number of candidates')
    ax.set_title('Number of cells supporting candidates')
    fig.tight_layout()
    if '/'.join(F1.split("/")[:-1]) == '':
        pltname = './{}_cells_supporting_candidate_indel.pdf'.format(PRES)
    else:
        pltname = '/'.join(F1.split("/")[:-1]) + "/{}_cells_supporting_candidate_indel.pdf".format(PRES)
    plt.savefig(pltname)


def get_readcounts(fin, loc):
    sample_dict = {}
    with open(fin, 'r') as fin:
        for line in fin:
            if not '+' in line:
                if not '-' in line:
                    continue
            line = line.strip().split()
            if line[0] + '_' + line[1] in loc:
                sample_dict[line[0] + '_' + line[1]] = {line[i]: line[i+1] for i in range(9, len(line), 2)}
                sample_dict[line[0] + '_' + line[1]]['RC'] = line[3]
    return sample_dict

def selectgoodcandidate(sample, candidate_readcount, snps_readcount, SAMPLES, CWD):
    from collections import defaultdict, deque
    import matplotlib.pyplot as plt
    import os

    sample_readcount = defaultdict(deque)
    for pos, v in candidate_readcount[sample].items():
        for s in SAMPLES:
            try:
                c = [c for k, c in snps_readcount[s][pos].items() if k==v[1]]
                if len(c) > 1: print('ERROR')
                if len(c) == 1: sample_readcount[s].append(int(c[0]))
                else: sample_readcount[s].append(0)
            except KeyError as e:
                sample_readcount[s].append(0)

    # Filter positions where the 'other' samples have as many alt reads as the focal samples
    bad_pos = deque()
    for i in range(len(sample_readcount[sample])):
        s_rc = sample_readcount[sample][i]
        s_all_rc = sum([sample_readcount[s][i] for s in SAMPLES])
        if s_all_rc >= 2*s_rc:
            bad_pos.append(i)
    bad_pos = set(bad_pos)
    # Select positions where the number of alt reads in other samples is less than
    # or equal to 20% of the alt reads count in focal sample
    good_pos = list(set(range(len(sample_readcount[sample]))) - bad_pos)
    high_conf = deque()
    for i in good_pos:
        s_rc = sample_readcount[sample][i]
        s_all_rc = sum([sample_readcount[s][i] for s in SAMPLES])
        if s_all_rc <= 1.2*s_rc:
            high_conf.append(i)

    # Get alt allele-freq at candidate positions in all samples
    sample_alfreq = defaultdict(deque)
    for pos, v in candidate_readcount[sample].items():
        for s in SAMPLES:
            try:
                c = [c for k, c in snps_readcount[s][pos].items() if k==v[1]]
                if len(c) > 1: print('ERROR')
                if len(c) == 1: sample_alfreq[s].append(int(c[0])/int(snps_readcount[s][pos]['RC']))
                else: sample_alfreq[s].append(0)
            except KeyError as e:
                sample_alfreq[s].append(0)
            except ZeroDivisionError as e:
                sample_alfreq[s].append(1)
    # Filter positions where the cumulative alt allele freq is more than or equal
    # to the focal sample
    bad_pos_af = deque()
    for i in range(len(sample_alfreq[sample])):
        s_rc = sample_alfreq[sample][i]
        s_all_rc = sum([sample_alfreq[s][i] for s in SAMPLES])
        if s_all_rc >= 2*s_rc:
            bad_pos_af.append(i)
    bad_pos_af = set(bad_pos_af)
    # Select positions where the alt allele-freq in other samples is less than or
    # equal to the alt allele freq in the focal sample
    good_pos_af = list(set(range(len(sample_alfreq[sample]))) - bad_pos_af)
    high_conf_af = deque()
    for i in good_pos_af:
        s_rc = sample_alfreq[sample][i]
        s_all_rc = (sum([sample_alfreq[s][i] for s in SAMPLES]) - s_rc)/3
        if s_all_rc <= 0.2*s_rc:
            high_conf_af.append(i)

    all_bad  = set(list(bad_pos) + list(bad_pos_af))
    all_high = set([i for i in high_conf if i in high_conf_af])
    all_mid  = set(range(len(sample_alfreq[sample]))) - set(list(all_bad) + list(all_high))

    with open('{}{}/{}_indel_uncertain_candidates.txt'.format(CWD, sample, sample), 'w') as fout:
        keys = list(candidate_readcount[sample].keys())
        for i in sorted(list(all_mid)):
            chr, p = keys[i].rsplit('_', 1)
            d = [str(c) for c in candidate_readcount[sample][keys[i]]]
            fout.write('\t'.join([chr, str(int(p)-1), p] + d) + '\n')
    with open('{}{}/{}_indel_good_candidates.txt'.format(CWD, sample, sample), 'w') as fout:
        keys = list(candidate_readcount[sample].keys())
        for i in sorted(list(all_high)):
            chr, p = keys[i].rsplit('_', 1)
            d = [str(c) for c in candidate_readcount[sample][keys[i]]]
            fout.write('\t'.join([chr, str(int(p)-1), p] + d) + '\n')

    positions_set = [set(range(len(sample_readcount[sample]))), all_bad, all_mid, all_high]
    pcnt = 1
    fig = plt.figure(figsize=[10, 10])
    marker = ['+', 'x', '.']
    plt.rc('font', size=8)
    for i in range(4):
        positions = positions_set[i]
        ax = fig.add_subplot(4, 2, pcnt)
        pcnt += 1
        m = 0
        for s in SAMPLES:
            if s == sample: continue
            ax.scatter([sample_readcount[sample][i] for i in positions], [sample_readcount[s][i] for i in positions], label=s, marker=marker[m])
            m += 1
        ax.legend(loc='upper left')
        ax.set_xlabel('{}'.format(sample))
        ax.set_ylabel('Other samples')
        if i == 0: title, lim = ('Alt allele readcount: All Candidates (#{})'.format(len(positions)), 250)
        elif i == 1: title, lim = ('Alt allele readcount: Bad Candidates (#{})'.format(len(positions)), 250)
        elif i == 2: title, lim = ('Alt allele readcount: Uncertain Candidates (#{})'.format(len(positions)), 250)
        elif i == 3: title, lim = ('Alt allele readcount: Good Candidates (#{})'.format(len(positions)), 250)
        ax.set_title(title)
        ax.set_ylim([-1, lim + 1])
        ax.set_xlim([-1, lim + 1])

        ax = fig.add_subplot(4, 2, pcnt)
        pcnt += 1
        m = 0
        for s in SAMPLES:
            if s == sample: continue
            ax.scatter([sample_alfreq[sample][i] for i in positions], [sample_alfreq[s][i] for i in positions], label=s, marker=marker[m])
            m += 1
        ax.legend(loc='upper left')
        ax.set_xlabel('{}'.format(sample))
        ax.set_ylabel('Other samples')
        if i == 0: title, lim = ('Alt allele frequency: All Candidates (#{})'.format(len(positions)), 1.05)
        elif i == 1: title, lim = ('Alt allele frequency: Bad Candidates (#{})'.format(len(positions)), 1.05)
        elif i == 2: title, lim = ('Alt allele frequency: Uncertain Candidates (#{})'.format(len(positions)), 1.05)
        elif i == 3: title, lim = ('Alt allele frequency: Good Candidates (#{})'.format(len(positions)), 1.05)
        ax.set_title(title)
        ax.set_ylim([-0.05, lim])
        ax.set_xlim([-0.05, lim])
    plt.tight_layout()
    plt.savefig('{}{}/{}_indel_alt_alleles.png'.format(CWD, sample, sample))
    plt.close()

def filterbackground(args):
    import os
    from collections import defaultdict, OrderedDict, deque
    from multiprocessing import Pool
    from functools import partial

    CWD = os.getcwd() + os.sep
    SAMPLES = ("MUT_11_1", "MUT_15", "WT_1", "WT_19")
    CORES = args.cores

    candidate_readcount = OrderedDict()
    for sample in SAMPLES:
        # print(sample)
        sample_dict = defaultdict()
        with open('{}/{}/multi_cell_{}_bt2_indel_candidate.sorted.filtered.bed'.format(CWD, sample, sample), 'r') as fin:
            for line in fin:
                line = line.strip().split()
                sample_dict[line[0] + '_' + line[2]] = (line[3], line[4], int(line[5]), float(line[6]))
        candidate_readcount[sample] = sample_dict
    candidate_readcount = OrderedDict(candidate_readcount)

    loc = set([k for sample in SAMPLES for k in candidate_readcount[sample].keys()])
    fins = ['{}{}/bam_read_counts_b30_q10.bt2.txt'.format(CWD, sample) for sample in SAMPLES]
    with Pool(processes=CORES) as pool:
        snps_readcount = pool.map(partial(get_readcounts, loc=loc), fins)
    snps_readcount = dict(zip(SAMPLES, snps_readcount))

    with Pool(processes=CORES) as pool:
        # selectgoodcandidate(sample, candidate_readcount, snps_readcount, CWD)
        pool.map(partial(selectgoodcandidate, candidate_readcount=candidate_readcount, snps_readcount=snps_readcount, SAMPLES=SAMPLES, CWD=CWD), SAMPLES)


if __name__=='__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()
    candidates = subparsers.add_parser("candidates", help="Get candidate indels in all samples", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    mcell = subparsers.add_parser("multicell", help="Select candidate indels that are present in multiple cells", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    filterbg = subparsers.add_parser("filterbg", help="Filter out candidates that are present in background samples", formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    candidates.set_defaults(func=getcandidate)
    candidates.add_argument('f', help='path to bam_readcount output files', type=argparse.FileType('r'), nargs='+')
    candidates.add_argument('-n', dest='n', help='minimum number of non-reference reads for selecting good candidates', type=int, default=5)
    candidates.add_argument('-m', dest='m', help='minimum coverage at position', type=int, default=50)
    candidates.add_argument('-M', dest='M', help='maximum coverage at position', type=int, default=350)
    candidates.add_argument('-s', dest='s', help='prefix for output files', type=str, default='pos', nargs='+')
    candidates.add_argument('-N', dest='N', help='Number of candidates to output (-1 for all candidates)', type=int, default=-1)
    candidates.add_argument('--cores', dest='C', help='Number of CPU cores to use', type=int, default=1)


    mcell.set_defaults(func=multicell)
    mcell.add_argument('pos', help='List of SNP candidate positions', type=argparse.FileType('r'))
    mcell.add_argument('rc_file_list', help='File containing paths to read count output files', type=argparse.FileType('r'))
    mcell.add_argument('-m', dest='m', help='minimum number of cells required to support a good candidate', type=int, default=3)
    mcell.add_argument('-o', dest='o', help='prefix for the output file', default='multi_cell', type=str)


    filterbg.set_defaults(func=filterbackground)
    filterbg.add_argument('--cores', help='Number of CPU cores to use', type=int, default=1)


    args = parser.parse_args()
    # print(args)
    # import sys
    # sys.exit()
    args.func(args)