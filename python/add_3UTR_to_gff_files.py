import argparse


def add3UTR(anno, fout, S):
    strand = anno[0][6]
    mrna_indices = [i for i in range(len(anno)) if anno[i][2] == 'mRNA'] + [len(anno)]
    hasutr = []
    for i in range(len(mrna_indices) - 1):
        utr3 = False
        for j in range(mrna_indices[i], mrna_indices[i+1]):
            if anno[j][2] == 'three_prime_UTR':
                utr3 = True
        hasutr.append(utr3)

    # All transcripts already have 3' UTR then print the annotation and return
    if False not in hasutr:
        for line in anno:
            fout.write('\t'.join(line) + '\n')
        return

    # For all transcripts that require 3' UTR, update exons and add 3'UTR.
    # Also, update gene coordinates if required.
    utranno = deque()
    for i in range(len(mrna_indices) - 1):
        if hasutr[i] == False:
            exon_indices = [j for j in range(mrna_indices[i], mrna_indices[i+1]) if anno[j][2] == 'exon']
            if strand=='+':
                anno[mrna_indices[i]][4] = str(int(anno[mrna_indices[i]][4]) + S)
                end = anno[mrna_indices[i]][4]
                utr3data = anno[exon_indices[-1]].copy()
                utr3data[2] = 'three_prime_UTR'
                utr3data[3] = str(int(utr3data[4]) + 1)
                utr3data[4] = str(int(utr3data[4]) + S)
                utr3data[8].replace(';Parent', 'utr3p1;Parent')
                utranno.append(utr3data)
                anno[exon_indices[-1]][4] = str(int(anno[exon_indices[-1]][4]) + S)
                if end != anno[exon_indices[-1]][4]:
                    print('ERROR: End coordinates do not match {}'.format(anno[0][8]))
            if strand == '-':
                anno[mrna_indices[i]][3] = str(int(anno[mrna_indices[i]][3])-S if (int(anno[mrna_indices[i]][3])-S) > 0 else 1)
                end = anno[mrna_indices[i]][3]
                utr3data = anno[exon_indices[-1]].copy()
                utr3data[2] = 'three_prime_UTR'
                utr3data[4] = str(int(utr3data[3])-1)
                utr3data[3] = str(int(utr3data[3])-S if (int(utr3data[3])-S) > 0 else 1)
                utr3data[8].replace(';Parent', 'utr3p1;Parent' )
                utranno.append(utr3data)
                anno[exon_indices[-1]][3] = str(int(anno[exon_indices[-1]][3]) - S if (int(anno[exon_indices[-1]][3]) - S) > 0 else 1)
                if end != anno[exon_indices[-1]][3]:
                    print('ERROR: End coordinates do not match {}'.format(anno[0][8]))
        else:
            utranno.append([])

    mincoord = min([int(anno[i][3]) for i in mrna_indices[:-1]])
    maxcoord = min([int(anno[i][4]) for i in mrna_indices[:-1]])
    anno[0][3] = str(mincoord)
    anno[0][4] = str(maxcoord)

    for i in range(len(mrna_indices)-2, -1, -1):
        if hasutr[i] == False:
            anno.insert(mrna_indices[i+1], utranno[i])
    for line in anno:
        fout.write('\t'.join(line) + '\n')
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Add 3'UTR regions to the input gff file",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gff', help='Input GFF file', type=argparse.FileType('r'))
    parser.add_argument('out', help='Output GFF file name', type=argparse.FileType('w'))
    parser.add_argument('-s', help='size of the UTR', type=int, default=501)

    args = parser.parse_args()
    GFF = args.gff.name
    OUT = args.out.name
    S = args.s
    S = 501

    from collections import deque
    geneid = ''
    strand = ''
    with open(GFF, 'r') as f:
        with open(OUT, 'w') as fout:
            anno = deque()
            for line in f:
                line = line.strip().split()
                if line[2] == 'gene':
                    if len(anno) > 0:
                        add3UTR(anno, fout, S)
                    anno = deque()
                    anno.append(line)
                else:
                    anno.append(line)
            add3UTR(anno, fout, S)