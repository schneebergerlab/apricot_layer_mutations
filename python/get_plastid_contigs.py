import argparse

def getplastidcontig(args):
    import sys
    sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/')
    from myUsefulFunctions import mergeRanges
    from collections import defaultdict
    from numpy import array

    paf = args.paf.name
    out = 'stdout' if args.o is None else args.o.name
    T = args.t
    contig_size = defaultdict(int)
    contig_map = defaultdict(list)
    with open(paf, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            contig_size[line[0]] = int(line[1])
            contig_map[line[0]].append([int(line[2]), int(line[3])])
    contig_map = {k: mergeRanges(array(v)) for k, v in contig_map.items()}
    contig_map_len = {k: sum(v[:, 1:]-v[:, :-1])[0]/contig_size[k] for k,v in contig_map.items()}
    if out == 'stdout': fout = sys.stdout
    else: fout = open(out, 'w')
    for k, v in contig_map_len.items():
        if v > T:
            fout.write('{}\t{}\n'.format(k, v))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('paf', help='Input alignments in PAF file format', type=argparse.FileType('r'))
    parser.add_argument('-o', dest='o', help='Output to file (default: STDOUT)', type=argparse.FileType('w'))
    parser.add_argument('-t', dest='t', help='Minimum contig length mapping to assign contig as plastid', type=float, default=0.8)
    args = parser.parse_args()
    getplastidcontig(args)



