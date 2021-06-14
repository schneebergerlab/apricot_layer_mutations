#!/usr/bin/env python3
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser("From the syri.out file, select SNP positions that are present in syntenic regions and do not overlap with SRs",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', help='path to syri.out file', type=argparse.FileType('r'))
    parser.add_argument('out', help='prefix for the output file', type=argparse.FileType('w'))

    args = parser.parse_args()

    from collections import deque

    snps = deque()
    sr = deque()
    with open(args.input.name, 'r') as fin:
    # with open('syri.out', 'r') as fin:
        for line in fin:
            line = line.strip().split()
            if line[3] in ['n', 'N'] or line[4] in ['n', 'N']: continue
            if line[10] == 'SNP' and 'SYN' in line[9]:
                snps.append([line[0], int(line[1]), line[3], line[4], line[5], int(line[6])])
                continue
            else:
                if line[10] in ['INV', 'TRANS', 'INVTR', 'DUP', 'INVDP']:
                    sr.append([line[0], int(line[1]), int(line[2]), line[5], int(line[6]), int(line[7])])

    from intervaltree import Interval, IntervalTree
    rtree = IntervalTree(Interval(r[1], r[2], r[0]) for r in list(sr))
    qtree = IntervalTree(Interval(r[4], r[5], r[3]) for r in list(sr))

    goodsnps = deque()
    for snp in snps:
        rng = rtree[snp[1]]
        if len(rng) == 0 or snp[0] not in [b.data for b in rng]:
            rng = qtree[snp[5]]
            if len(rng) == 0 or snp[4] not in [b.data for b in rng]:
                goodsnps.append(snp)

    with open(args.out.name, 'w') as fout:
        for snp in goodsnps:
            fout.write(" ".join([snp[0], str(snp[1]), str(snp[1]), snp[2], snp[3]]) + "\n")




