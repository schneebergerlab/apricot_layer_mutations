# from Bio.Seq import reverse_complement
from multiprocessing import Pool
from get_primer_match import getmatch
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('primer', help='Primer sequence to check', type=str)
    parser.add_argument('maxmm', help='Maximum allowed mismatches', type=int, default=4)
    parser.add_argument('-o', help='Output file name', type=argparse.FileType('w'))
    parser.add_argument('-n', help='number of cores to use', default=1, type=int)

    args = parser.parse_args()
    PRIMER = args.primer
    M = args.maxmm
    O = args.o.name if args.o is not None else 'primer_match.txt'
    NCORES = args.n

    import sys
    sys.path.insert(0, '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/primer_match_find')
    from myUsefulFunctions import readfasta

    curg = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/manual_curated_assembly/cur/cur.genome.v1.fasta'
    orag = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/manual_curated_assembly/ora/ora.genome.v1.fasta'

    cur = readfasta(curg)
    ora = readfasta(orag)
    comparisons = [[PRIMER, s, M, c] for gen in [cur, ora] for c, s in gen.items()]
    with Pool(processes=NCORES) as pool:
    # with Pool(processes=3) as pool:
        matches = pool.starmap(getmatch, comparisons)

    with open(O, 'w') as fout:
        for match in matches:
            for m in match:
                fout.write('\t'.join(m) + '\n')
