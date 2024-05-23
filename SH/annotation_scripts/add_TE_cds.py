"""
Read a GFF file in which transposable element annatations have exons but not cds
and add CDS information to the file
"""
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser("Adds CDS info to the transposable element annotation in the GFF file.",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gff', help='Input gff file', type=argparse.FileType('r'))
    parser.add_argument('out', help='Output gff file', type=argparse.FileType('w'))

    args = parser.parse_args()

    with open(args.out.name, 'w') as fout:
        with open(args.gff.name, 'r') as fin:
            for line in fin:
                fout.write(line)
                if line[0] == '#': continue
                if line.strip().split()[2] == 'gene':
                    if 'trans' in line:
                        trans = True
                        count = 1
                    else: trans = False
                    continue
                if line.strip().split()[2] in ['mRNA', 'CDS']: continue
                if line.strip().split()[2] != 'exon':
                    sys.exit('Unknown annotations type in line {}'.format(line))
                if trans:
                    line = line.strip().split()
                    line[2] = 'CDS'
                    parent = line[8].split(';')[1].split('=')[1]
                    line[8] = 'ID={}.cds{};Parent={}'.format(parent, count, parent)
                    fout.write('\t'.join(line) + '\n')
                    count += 1



