import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser("Get all mutations present in a given cell",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('f1', help='List of mutations', type=argparse.FileType('r')) # in the format generated by /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/select_candidate_supported_by_cells.py
    parser.add_argument('f2', help='bam-readcount output for the cell', type=argparse.FileType('r'))
    parser.add_argument('bc', help='barcode of the cell', type=str)
    parser.add_argument('-n', dest='n', help='minimum alternate read count to select position', type=int, default=1)
    parser.add_argument('-o', dest='o', help='prefix for the output file', default='cell_var', type=str)

    args = parser.parse_args()

    BASE_DICT = {
        6: "A",
        7: "C",
        8: "G",
        9: "T"
    }
    F1      = args.f1.name
    F2      = args.f2.name
    N       = args.n
    PRES    = args.o






