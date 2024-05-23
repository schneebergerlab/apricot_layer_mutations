#!/usr/bin/env python3
import argparse

if __name__ == '__main__':
    """
    Get the CDS regions start and end points in BED format from the GFF file.
    
    
    """
    parser = argparse.ArgumentParser("In PASA output gff file, select single isoform for a gene",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('input', help='Input GFF file', type=argparse.FileType('r'))
    parser.add_argument('out', help='Output BED file Name', type=argparse.FileType('w'))
    parser.add_argument('sample', help='Sample Name', type=str)

    args = parser.parse_args()
    
    from collections import deque
    import sys
    with open(args.input.name, 'r') as fin:
        with open(args.out.name, 'w') as fout:
            mRNAid = ''
            geneID = ''
            cds = deque()
            for line in fin:
                line = line.strip().split()
                if line[2] == 'gene':
                    if geneID == '':
                        geneID = line[8].split(';')[0].split('=')[1]
                        chr = line[0]
                        continue
                    mRNAid = mRNAid.replace('mRNA', args.sample)
                    fout.write('\t'.join([chr, str(min(cds)), str(max(cds)), mRNAid]) + '\n')
                    geneID = line[8].split(';')[0].split('=')[1]
                    chr = line[0]
                    cds = deque()
                    continue
                if line[2] == 'mRNA':
                    parent = line[8].split(';')[1].split('=')[1]
                    if parent != geneID: sys.exit("Error: mRNA do not match gene {}".format(geneID))
                    mRNAid = line[8].split(';')[0].split('=')[1]
                    continue
                if line[2] != 'CDS': continue
                cds.append(int(line[3]))
                cds.append(int(line[4]))
            mRNAid = mRNAid.replace('mRNA', args.sample)
            fout.write('\t'.join([chr, str(min(cds)), str(max(cds)), mRNAid]) + '\n')

    