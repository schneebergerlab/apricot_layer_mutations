import argparse

if __name__ == '__main__':
    """
    For the candidate mutation positions, count the number of cells that have reads supporting the mutation.
    """

    parser = argparse.ArgumentParser("Get number of cells supporting mutations",
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('pos', help='List of SNP candidate positions', type=argparse.FileType('r'))
    parser.add_argument('pos_info', help='Read count information for the positions', type=argparse.FileType('r'))
    parser.add_argument('rc_file_list', help='File containing paths to read count output files', type=argparse.FileType('r'))
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
    F1 = args.pos.name
    F2 = args.pos_info.name
    F3 = args.rc_file_list.name
    N  = args.n
    M  = args.m
    PRES  = args.o



    import os
    import pandas as pd
    from collections import defaultdict, deque, Counter

    # if not os.path.isdir(D1):
        # sys.exit('Incorrect path to cell read count folder')

    f1 = pd.read_table(F1, header=None)
    f2 = pd.read_table(F2, header=None, sep=' ')
    f3 = pd.read_table(F3, header=None)


    # f1 = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_1/WT_1_only_SNPs_candidate.sorted.common.regions', header=None)
    f1.columns = ['chr', 'pos0', 'pos', 'vac', 'var']
    
   
    # f2 = pd.read_table('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_1/WT_1_only_SNPs_candidate.txt',
    #                    header=None, sep=' ')
    f2.sort_values([2, 3], inplace=True)
    f2.columns = ['vac', 'var', 'chr', 'pos', 'ref', 'rc', 'A', 'C', 'G', 'T', 'N', 'vac_rank', 'var_rank', 'rank']
    
    
    # f3 = f3.loc[f3['vac'] >= N]
    # f3.index = range(f3.shape[0])
    f3.columns = ['fin']

    # f1 = f1.merge(f3, on=['chr', 'pos'], how='inner')
    # f2 = f2.merge(f3, on=['chr', 'pos'], how='inner')

    alt = defaultdict(dict)
    # f2_alt = defaultdict(dict)
    for row in f2.itertuples(index=False):
        for i in range(6, 10):
            if BASE_DICT[i] == row.ref: continue
            if row[i] >= N:
                alt[row.chr][row.pos] = (row.ref, BASE_DICT[i])
                break

    # for row in f2.itertuples(index=False):
        # for i in range(6, 10):
            # if BASE_DICT[i] == row.ref: continue
            # if row[i] >= N:
                # f2_alt[row.chr][row.pos] = (row.ref, BASE_DICT[i])
                # break
                

    f1_alts = deque()
    f1_refs = deque()
    for row in f1.itertuples(index=False):
      f1_refs.append(alt[row.chr][row.pos][0])
      f1_alts.append(alt[row.chr][row.pos][1])
     
    f1['ref'] = f1_refs
    f1['alt'] = f1_alts
    # bad_candidates = [i for i in range(len(f1_alt)) if f1_alt[i] != f2_alt[i]]
    # bad_locs = set(f3.iloc[bad_candidates]['chr'] + '_' + f3.iloc[bad_candidates]['pos'].astype(str))
    # f3 = f3.drop(bad_candidates, axis=0)
    #
    # bad_locs = [row.chr + "_" + str(row.pos) for row in f3.itertuples(index=False) if row.alt == 'BAD']
    # f3 = f3.loc[f3.ref != 'BAD']
    f1_alt = {row.chr + '_' + str(row.pos): [row.alt, 0] for row in f1.itertuples(index=False)}

    BASE_DICT2 = {
        "A": 4,
        "C": 5,
        "G": 6,
        "T": 7
    }
    # bcs = sorted(os.listdir(D1))
    for bc in f3.fin:
        # with open(D1 + os.sep + bc + os.sep + bc + '_candidate_readcount.txt') as fin:
        with open(bc) as fin:
            for line in fin:
                line = line.strip().split()
                # if line[0] + '_' + line[1] in bad_locs: continue
                if int(line[BASE_DICT2[f1_alt[line[0] + '_' + line[1]][0]]]) >= 1:
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
        pltname = './{}_cells_supporting_candidate.png'.format(PRES)
    else:
        pltname = '/'.join(F1.split("/")[:-1]) + "/{}_cells_supporting_candidate.png".format(PRES)
    plt.savefig(pltname)
