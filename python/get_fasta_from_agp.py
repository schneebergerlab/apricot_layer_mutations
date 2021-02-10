"""
Temporary which accepts manually edited file as input and generates fasta genome from it.
Contigs are joined in the order in which they appear in the agp file
POTENTIALLY: can be generated as a general purpose method.
"""

import argparse


def getfasta(args):
    from collections import OrderedDict, defaultdict
    import sys
    sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/')
    from myUsefulFunctions import readfasta, writefasta
    from Bio.Seq import reverse_complement

    AGP = args.agp.name
    OUT = 'stdout' if args.o is None else args.o
    FIN = args.fasta.name
    Ns = ''.join(['N']*args.N)

    agpdata = defaultdict(OrderedDict)
    with open(AGP, 'r') as fin:
        for line in fin:
            if line[0] == '#': continue
            line = line.strip().split()
            if line[4] != 'W': continue
            agpdata[line[0]][line[5]] = line[8]

    fasta = readfasta(FIN)

    outfasta = OrderedDict()
    for k, v in agpdata.items():
        out_seq = [fasta[cid] if strand == '+' else reverse_complement(fasta[cid]) for cid, strand in v.items()]
        outfasta[k] = Ns.join(out_seq)
    writefasta(outfasta, OUT)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('agp', help='Input .agp file to create fasta', type=argparse.FileType('r'))
    parser.add_argument('fasta', help='Input fasta file containing contigs', type=argparse.FileType('r'))
    parser.add_argument('-o', dest='o', help='Output fasta file name', type=str, default='stdout')
    parser.add_argument('-N', dest='N', help='Number of Ns to insert between two contigs', type=int, default=100)

    args=parser.parse_args()
    getfasta(args)