#!/usr/bin/env python
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("snps", help="SNP positions", type=argparse.FileType('r'))
    parser.add_argument("rc", help="Readcount at SNP positions", type=argparse.FileType('r'))
    parser.add_argument("out", help="Output file name", type=argparse.FileType('w'))

    args = parser.parse_args()

    SNPFIN = args.snps.name
    RCFIN = args.rc.name
    OUTFIN = args.out.name
    BASEDICT = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
    from collections import defaultdict, deque



    snps = defaultdict(dict)
    # SNPFIN = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/strict_syn_snp.txt'
    with open(SNPFIN, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            snps[line[0]][line[1]] = deque([[line[3], line[4]]])

    # RCFIN = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/strict_syn_snp_readcount.txt'
    with open(RCFIN, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            r = snps[line[0]][line[1]][0][0]
            q = snps[line[0]][line[1]][0][1]
            snps[line[0]][line[1]].append([int(line[BASEDICT[r]]), int(line[BASEDICT[q]])])

    # OUTFIN='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/TMP_snps_readcount_stats.txt'
    with open(OUTFIN, 'w') as fout:
        for chr, poss in snps.items():
            for pos, v in poss.items():
                # print(chr, pos, v, t, maf)
                try:
                    t = sum(v[1])
                    try:
                        maf = v[1][1]/t
                    except ZeroDivisionError:
                        maf = 0
                    fout.write('\t'.join([chr, pos, v[0][0], str(v[1][0]), v[0][1], str(v[1][1]), str(t), str(round(maf, 3))]) + '\n')
                except IndexError:
                    t = 0
                    maf = 0
                    fout.write('\t'.join([chr, pos, v[0][0], '0', v[0][1], '0', '0', '0']) + '\n')
                # print(chr, pos, v, t, maf)
