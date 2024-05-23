#!/usr/bin/env python3
"""
Read the input GFF file and select only one (first) isoform for each gene.
This script is written so that output of PASA which contains multiple isoforms
can be used easily with the down-stream manual curation.
"""

import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser("In PASA output gff file, select single isoform for a gene",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('input', help='Input GFF file', type=argparse.FileType('r'))
    parser.add_argument('out', help='Output GFF file', type=argparse.FileType('w'))    
    
    args = parser.parse_args()
    # print(args.in)
    import sys
    with open(args.input.name, 'r') as fin:
        with open(args.out.name, 'w') as fout:
            mRNAid = ''
            geneID = ''
            for line in fin:
                line = line.strip().split()
                if line[2] == 'gene':
                    geneID = line[8].split(';')[0].split('=')[1]
                    fout.write("\t".join(line) + '\n')
                    continue
                if line[2] == 'mRNA':
                    parent = line[8].split(';')[1].split('=')[1]
                    if parent != geneID: sys.exit("Error: mRNA do not match gene {}".format(geneID))
                    mRNAid = line[8].split(';')[0].split('=')[1]
                    if mRNAid[-2:] != '.1':
                        mRNAid = ''
                        continue
                    else:
                        fout.write("\t".join(line) + '\n')
                        continue
                if mRNAid == '': continue
                fout.write("\t".join(line) + '\n')



