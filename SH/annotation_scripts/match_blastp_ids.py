#!/usr/bin/env python3
import argparse

if __name__ == '__main__':
    """
    Revert back gene ids for the blastp output from orthofinder
    """
    parser = argparse.ArgumentParser("Map gene IDs in the blastp output file to the original gene ids",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('blastp_out', help='Output blastp file from orthofinder', type=argparse.FileType('r'))
    parser.add_argument('sequence_id', help='SequenceIDs.txt from orthofinder', type=argparse.FileType('r'))
    parser.add_argument('out', help='Output blastp file', type=argparse.FileType('w'))

    args = parser.parse_args()
    
    from collections import deque
    import sys
    sequence_ids = {}
    with open(args.sequence_id.name, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            sequence_ids[line[0][:-1]] = line[1]
        
    with open(args.blastp_out.name, 'r') as fin:
        with open(args.out.name, 'w') as fout:
            for line in fin:
                line = line.strip().split()
                line[0] = sequence_ids[line[0]]
                line[1] = sequence_ids[line[1]]
                fout.write('\t'.join(line) + '\n')
    