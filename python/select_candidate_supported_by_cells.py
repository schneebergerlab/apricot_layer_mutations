import argparse

if __name__ == '__main__':
    """
    For the candidate mutation positions, count the number of cells that have reads supporting the mutation.
    """

    parser = argparse.ArgumentParser("Get number of cells supporting mutations",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('f1', help='path to minimap2 candidate list', type=argparse.FileType('r'))
    parser.add_argument('f2', help='path to bowtie2 candidate list', type=argparse.FileType('r'))
    parser.add_argument('f3', help='path to common candidate bed file', type=argparse.FileType('r'))
    parser.add_argument('d1', help='path to folder containing cells read count', type=str)
    parser.add_argument('-n', dest='n', help='minimum read count used to select candidates', type=int, default=4)
    parser.add_argument('-m', dest='m', help='minimum number of cells required to support a good candidate', type=int, default=3)
    parser.add_argument('-o', dest='o', help='prefix for the output file', default='multi_cell', type=str)

    args = parser.parse_args()
    BASE_DICT = {
        6: "A",
        7: "C",
        8: "G",
        9: "T"
    }
    F1 = args.f1.name
    F2 = args.f2.name
    F3 = args.f3.name
    D1 = args.d1
    N  = args.n
    M  = args.m
    PRES  = args.o


    # D1 = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/MUT_11_1/cells_readcount'

    import os
    if not os.path.isdir(args.d1): sys.exit('Incorrect path to cell read count folder')

    import pandas as pd
    f1 = pd.read_table(F1, header=None, sep=' ')
    f2 = pd.read_table(F2, header=None, sep=' ')
    f3 = pd.read_table(F3, header=None)


    # f1 = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_19/WT_19_only_SNPs_candidate.txt',
    #                   header=None, sep=' ')
    f1.sort_values([2, 3], inplace=True)
    f1.columns = ['vac', 'var', 'chr', 'pos', 'ref', 'rc', 'A', 'C', 'G', 'T', 'N', 'vac_rank', 'var_rank', 'rank']

    # f2 = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_19/WT_19_bt2_only_SNPs_candidate.txt',
    #                   header=None, sep=' ')
    f2.sort_values([2, 3], inplace=True)
    f2.columns = ['vac', 'var', 'chr', 'pos', 'ref', 'rc', 'A', 'C', 'G', 'T', 'N', 'vac_rank', 'var_rank', 'rank']

    # f3 = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_19/WT_19_only_SNPs_candidate.sorted.common.regions', header=None)
    f3.columns = ['chr', 'pos0', 'pos', 'vac', 'var']
    f3 = f3.loc[f3['vac'] >= N]
    f3.index = range(f3.shape[0])

    f1 = f1.merge(f3, on=['chr', 'pos'], how='inner')
    f2 = f2.merge(f3, on=['chr', 'pos'], how='inner')

    from collections import deque
    f1_alt = deque()
    f2_alt = deque()
    count = 0
    for row in f1.itertuples(index=False):
        for i in range(6, 10):
            if BASE_DICT[i] == row.ref: continue
            if row[i] >= N:
                count += 1
                f1_alt.append(BASE_DICT[i])
                break

    for row in f2.itertuples(index=False):
        for i in range(6, 10):
            if BASE_DICT[i] == row.ref: continue
            if row[i] >= N:
                f2_alt.append(BASE_DICT[i])
                break

    from collections import Counter
    assert(all(f1.ref == f2.ref))
    f3['ref'] = f1.ref
    f3['alt'] = list(f1_alt)
    bad_candidates = [i for i in range(len(f1_alt)) if f1_alt[i] != f2_alt[i]]
    bad_locs = set(f3.iloc[bad_candidates]['chr'] + '_' + f3.iloc[bad_candidates]['pos'].astype(str))
    f3 = f3.drop(bad_candidates, axis=0)

    f3_alt = {row.chr + '_' + str(row.pos): [row.alt, 0] for row in f3.itertuples(index=False)}

    BASE_DICT2 = {
        "A": 4,
        "C": 5,
        "G": 6,
        "T": 7
    }
    bcs = sorted(os.listdir(D1))
    for bc in bcs:
        # print(bc)
        with open(D1 + os.sep + bc + os.sep + bc + '_candidate_readcount.txt') as fin:
            for line in fin:
                line = line.strip().split()
                if line[0] + '_' + line[1] in bad_locs: continue
                if int(line[BASE_DICT2[f3_alt[line[0] + '_' + line[1]][0]]]) >= 1:
                    f3_alt[line[0] + '_' + line[1]][1] += 1

    f3_loci = list(f3['chr'].astype(str) + '_' + f3['pos'].astype(str))
    f3_g = [True if f3_alt[i][1] >= 3 else False for i in f3_loci]

    f3_g_can = f3.loc[f3_g].copy()

    if '/'.join(F3.split("/")[:-1]) == '':
        fout = './' + "/{}".format(PRES) + '_' + F3.split("/")[-1]
    else:
        fout = '/'.join(F3.split("/")[:-1]) + "/{}".format(PRES) + '_' + F3.split("/")[-1]
    f3_g_can.to_csv(path_or_buf=fout,
              sep='\t',
              header=False,
              index=False)


    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=[5, 3])
    ax = fig.add_subplot(1, 1, 1)
    ax.bar(x= list(Counter([i[1] for i in f3_alt.values()]).keys()),
           height=list(Counter([i[1] for i in f3_alt.values()]).values()))
    ax.set_yscale('log')
    ax.set_xlabel('Number of cells')
    ax.set_ylabel('Number of candidates')
    ax.set_title('Number of cells supporting candidates')
    fig.tight_layout()

    if '/'.join(F3.split("/")[:-1]) == '':
        pltname = './cells_supporting_candidate.png'
    else:
        pltname = '/'.join(F3.split("/")[:-1]) + "/cells_supporting_candidate.png"
    plt.savefig(pltname)