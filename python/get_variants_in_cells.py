#!/usr/bin/env python
import argparse

BASE_DICT = {
    "A": 4,
    "C": 5,
    "G": 6,
    "T": 7
}

def get_rc(f, muts):
    with open(f, 'r') as fin:
        keys = list(muts.keys())
        keys2 = set(muts.keys())
        ref_rc = ['0']*len(keys)
        alt_rc = ['0']*len(keys)
        for line in fin:
            line = line.strip().split()
            if line[0] + '_' + line[1] not in keys2: continue
            if int(line[3]) == 0: continue
            i = keys.index(line[0] + '_' + line[1])
            ref_rc[i] = line[BASE_DICT[muts[line[0] + '_' + line[1]][0]]]
            alt_rc[i] = line[BASE_DICT[muts[line[0] + '_' + line[1]][1]]]
    return (ref_rc, alt_rc)



if __name__ == '__main__':
    parser = argparse.ArgumentParser("Get allele depth of mutations in read-count data",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('mutations', help='List of mutations', type=argparse.FileType('r')) # BED-like file. For read position only the third column is used.
    parser.add_argument('-rc', help='bam-readcount output for the cell', type=argparse.FileType('r'))
    parser.add_argument('-frc', help='File containing paths to bam-readcount files. One per line', type=argparse.FileType('r'))
    # parser.add_argument('-n', dest='n', help='minimum alternate read count to select position', type=int, default=1)
    parser.add_argument('-o', dest='o', help='prefix for the output file', default='mut_rc', type=str)

    args = parser.parse_args()

    import sys


    F1      = args.mutations.name
    # N       = args.n
    PRES    = args.o

    if args.rc is not None and args.frc is not None:
        raise ValueError('-rc and -frc both cannot be provided. Using only -frc')
    if args.rc is None and args.frc is None:
        raise ValueError('Need bam-readcount path. Both -rc and -frc are not defined. Exiting.')
        sys.exit()
    if args.rc is not None: F2      = args.rc.name
    else : F2      = args.frc.name

    from collections import OrderedDict
    muts = OrderedDict()
    with open(F1, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            muts[line[0] + '_' + line[2]] = (line[3], line[4], line[5], line[6])

    keys = list(muts.keys())
    if args.rc is not None:
        ref_rc, alt_rc = get_rc(F2, muts)
        with open(PRES + ".txt", 'w') as fout:
            fout.write('\t'.join(["chr", "start", "end", "ref", "alt", "depth", "freq"]) + '\t' + '\t'.join([F2+"_ref", F2+"_alt"]) + '\n')
            for i in range(len(keys)):
                k = keys[i].split('_')
                v = muts[keys[i]]
                fout.write('\t'.join([k[0], k[1], k[1], v[0], v[1], v[2], v[3], ref_rc[i], alt_rc[i]]) + '\n')
    else:
        t = ['ref', 'alt']
        with open(F2, 'r') as fin:
            rcs = OrderedDict()
            for line in fin:
                line = line.strip().split()
                rcs[line[1]] = get_rc(line[0], muts)
                print(line[1])
            samples = list(rcs.keys())

            with open(PRES + ".txt", 'w') as fout:
                fout.write('\t'.join(["chr", "start", "end", "ref", "alt", "depth", "freq"]) + '\t' + '\t'.join([s+';'+t[j] for s in samples for j in range(2)]) + '\n')
                for i in range(len(keys)):
                    k = keys[i].split('_')
                    v = muts[keys[i]]
                    fout.write('\t'.join([k[0], k[1], k[1], v[0], v[1], v[2], v[3]] + [rcs[s][j][i] for s in samples for j in range(2)]) + '\n')











