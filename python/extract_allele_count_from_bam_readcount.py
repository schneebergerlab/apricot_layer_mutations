#!/usr/bin/env python3

import argparse

if __name__ == '__main__':
    '''
    SNP Info Columns:   Chr Pos Pos Ref Alt
    Readcount Columns:  Chr Pos Ref Count A_cnt C_cnt G_cnt T_cnt  
    '''
    parser = argparse.ArgumentParser('Get allele counts at given marker position from the bam-readcount output file')
    parser.add_argument(dest='rc', help='Path to read count file', type=argparse.FileType('r'))
    parser.add_argument(dest='snp', help='Path to high-quality/corrected markers list file', type=argparse.FileType('r'))
    parser.add_argument('o', help='Output file name', type=argparse.FileType('w'))

    args = parser.parse_args()
    # print(args)
    # print(args.hqfin.name)
    #
    import pandas as pd
    from collections import deque, defaultdict

    base2pos = {'A': 4,
                'C': 5,
                'G': 6,
                'T': 7
                }

    hqmarkers = dict()
    badpos = set()
    with open(args.snp.name, 'r') as fin:
    # with open('../../strict_syn_snp.selected.txt', 'r') as fin:
        for line in fin:
            line = line.strip().split()
            if line[0]+line[1] in badpos: continue
            if line[0]+line[1] in hqmarkers:
                if hqmarkers[line[0]+line[1]] != (line[3], line[4]):
                    hqmarkers.pop(line[0]+line[1])
                    badpos.add(line[0]+line[1])
            else:
                hqmarkers[line[0]+line[1]] = (line[3], line[4])

    keys = set(hqmarkers.keys())
    rc = deque()
    with open(args.rc.name, 'r') as fin:
    # with open('AAACCTGAGCCCAACC_read_counts_b30_q40.bt2.txt', 'r') as fin:
        for line in fin:
            line = line.strip().split()
            if line[0]+line[1] in keys:
                if line[3] == '0': continue
                rc.append([line[0], line[1], hqmarkers[line[0]+line[1]][0], line[base2pos[hqmarkers[line[0]+line[1]][0]]],  hqmarkers[line[0]+line[1]][1], line[base2pos[hqmarkers[line[0]+line[1]][1]]]])
                keys.remove(line[0]+line[1])

    rcdf = pd.DataFrame(rc)
    rcdf[1] = rcdf[1].astype(int)
    rcdf.sort_values([0, 1], inplace=True)
    rcdf.to_csv(args.o.name, header=False, index=False, sep='\t')
    # rcdf.to_csv('TMP.txt', header=False, index=False, sep='\t')

