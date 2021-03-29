import argparse

def get_pos(fin):
    with open(fin, 'r') as f:
        out = defaultdict(deque)
        count = 0
        for line in f:
            count += 1
            line = line.strip().split()
            out[line[0]].append(line[1])
    return out



def write_pos_only_snps(fin, fout, pos, CAN_N, RC_MIN, RC_MAX, CAN_CNT):
    '''
    Select all positions that one alternate substitution allele (SNPs).
    '''
    base_dict = {4: 'A',
                 5: 'C',
                 6: 'G',
                 7: 'T'
                 }
    
    with open(fin, 'r') as f:
        df = deque()
        for line in f:
            line = line.strip().split()
            if int(line[3]) < RC_MIN: continue
            if int(line[3]) > RC_MAX: continue
            hf = 0
            # for i in [4, 5, 6, 7, 10, 12, 14]:
            for i in [4, 5, 6, 7]:
                try:
                    if int(line[i]) >= CAN_N:
                        hf += 1
                except IndexError as e:
                    break
            if hf > 2:
                continue
            try:
                if line[1] in pos[line[0]]:
                    ## Get alternate allele readcount and alternate allele ratio
                    rc = 0
                    ar = 0
                    alt = ''
                    for i in [4, 5, 6, 7]:
                        try:
                            if int(line[i]) >= CAN_N:
                                if base_dict[i] != line[2]:
                                    rc = int(line[i])
                                    ar = int(line[i])/int(line[3])
                                    alt = base_dict[i]
                        except IndexError as e:
                            pass
                    if rc == 0 or ar == 0 or alt == '':
                        continue
                    df.append([line[0], line[1], line[1], line[2], alt, rc, ar])
            except KeyError as e:
                pass

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


if __name__=='__main__':
    """
    from the filtered bam-readcount output, use different thresholds to find candidate mutations
        ## Conditions for filtering:
        # Read count should be less than a given cutoff
        # The alternate allele should have at least three reads supporting it
        # Only one alternate should have more than 2 reads
        # Same position is not variant in other branches

        # Secondary task
        ## Perform on both (MM2 and BT2) alignments: and then select positions which are selected by both methods
    """

    parser = argparse.ArgumentParser("Get candidate somatic mutation positions",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('f', help='path to bam_readcount output files', type=argparse.FileType('r'), nargs='+')
    parser.add_argument('-n', dest='n', help='minimum number of non-reference reads for selecting good candidates', type=int, default=4)
    parser.add_argument('-m', dest='m', help='minimum coverage at position', type=int, default=50)
    parser.add_argument('-M', dest='M', help='maximum coverage at position', type=int, default=350)
    parser.add_argument('-s', dest='s', help='prefix for output files', type=str, default='pos', nargs='+')
    parser.add_argument('-N', dest='N', help='Number of candidates to output (-1 for all candidates)', type=int, default=-1)

    args = parser.parse_args()

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

    from collections import deque, defaultdict
    import pandas as pd
    from matplotlib import pyplot as plt
    from upsetplot import UpSet


    pos_list = deque()
    for f in FINS:
        pos_list.append(get_pos(f))

    for i in range(N):
        for k in pos_list[i].keys():
            pos_list[i][k] = set(pos_list[i][k])

    for i in range(N):
        # For each sample select and write only SNP positions
        if '/'.join(FINS[i].split("/")[:-1]) == '':
            fout = './' + "/{}".format(PRES[i]) + 'SNPs_candidate'
        else:
            fout = '/'.join(FINS[i].split("/")[:-1]) + "/{}".format(PRES[i]) + 'SNPs_candidate'
        write_pos_only_snps(FINS[i], fout, pos_list[i], CAN_N, RC_MIN, RC_MAX, CAN_CNT)

    ## Reading Bed files generated above to perform filtersing
    import pybedtools
    beds = deque()
    for i in range(N):
        if '/'.join(FINS[i].split("/")[:-1]) == '':
            fout2 = './' + "/{}".format(PRES[i]) + 'SNPs_candidate'
        else:
            fout2 = '/'.join(FINS[i].split("/")[:-1]) + "/{}".format(PRES[i]) + 'SNPs_candidate'
        beds.append(pybedtools.BedTool(fout2 + '.sorted.bed'))


    ## Generate UpSet plot for the mutations: Uncomment to use this sections

    # beds_df = deque()
    # for i in range(N):
        # pos = deque()
        # for b in beds[i]:
            # # if int(list(b)[5]) < 10: continue
            # pos.append('_'.join(list(b)[:2]))
        # pos = pd.DataFrame(list(pos))

        # # pos[1] = pos[1].astype('int')
        # beds_df.append(pos)

    # for i in range(N):
        # beds_df[i].sort_values([0], inplace=True)

    # merge_df = beds_df[0]
    # for i in range(1, N):
        # merge_df = merge_df.merge(beds_df[i], 'outer', on = [0])

    # merge_df.sort_values([0], inplace=True)
    # merge_df.index = range(merge_df.shape[0])

    # for i in range(N):
        # merge_df[args.s[i]] = merge_df[0].isin(beds_df[i][0])
    # merge_df = merge_df.set_index(args.s)
    # upset = UpSet(merge_df, subset_size='count', intersection_plot_elements=3, element_size=40)
    # fig = plt.figure(figsize=[8, 8])
    # upset.plot(fig=fig)
    # plt.yscale('log')
     # # plt.show()
    # plt.savefig('_'.join(args.s) + ".pdf")



    """
    Earlier, I was first selecting positions unique in all samples and then was
    finding for candidate positions. That was too stringent filtering as noise
    in other samples could also filter out signal from the focal sample.

    Now, I first select candidate SNPs in each sample, and then select sample
    specific SNPs by removing candidate SNPs in other samples.
    """
    ## Select positions that are mutated in all samples
    beds_reg = {}
    for i in range(N):
        # count = 0
        for b in beds[i]:
            # count += 1
            # if count == 101: break
            try:
                beds_reg['_'.join(list(b)[:5])] += 1
            except KeyError  as e:
                beds_reg['_'.join(list(b)[:5])] = 1
            # break
    bad_pos = set([k for k, v in beds_reg.items() if v == N])


    # inallbeds = beds[0] + beds[1] + beds[2] + beds[3]
    # beds_uni = beds.copy()
    # for i in range(N):
        # # rem_bed = set(range(N)) - set([i])
        # # for j in rem_bed:
        # #     beds_uni[i] = beds_uni[i] - beds[j]
        # beds_uni[i] = beds_uni[i] - inallbeds


    ## Output positions that are not mutated in all samples
    for i in range(N):
        if '/'.join(FINS[i].split("/")[:-1]) == '':
            fout = './' + "/{}".format(PRES[i]) + 'filtered_SNPs_candidate.sorted.bed'
        else:
            fout = '/'.join(FINS[i].split("/")[:-1]) + "/{}".format(PRES[i]) + 'filtered_SNPs_candidate.sorted.bed'
        with open(fout, 'w') as f:
            for b in beds[i]:
                if '_'.join(list(b)[:5]) not in bad_pos:
                    f.write('\t'.join(list(b)) + '\n')
        # beds_uni[i].saveas(fout)
