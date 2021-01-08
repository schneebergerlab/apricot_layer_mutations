import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser("Add 3'UTR regions to the input gff file",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gff', help='Input GFF file', type=argparse.FileType('r'))
    parser.add_argument('-s', help='size of the UTR', type=int, default=500)
    parser.add_argument('-o', help='output file name', type=str, default='out.gff')

    args = parser.parse_args()
    S = args.s
    S = 500

    geneid = ''
    strand = ''
    end    = ''
    with open(args.gff.name, 'r') as f:
    with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/Cur.v1.1.protein.coding.genes.gff', 'r') as f:
        with open('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/Cur.v1.1.protein.coding.genes.manually_edited.gff', 'w') as fout:
            for line in f:
                line = line.strip().split()
                # break
                # print(line)
                if line[2] == 'gene':
                    strand = line[6]

                    if strand == '+':
                        end = int(line[4])
                        line[4] = str(int(line[4]) + S)
                    elif strand == '-':
                        end = int(line[3])
                        line[3] = str(max(0, int(line[3]) - S))
                    else: raise ValueError("Incorrect value of strand.\n" + '\t'.join(line))

                    geneid = line[8].split(';')[0].split('=')[1]


                else:
                    if line[8].split(';')[1].split('=')[1].split('.')[0] != geneid: raise ValueError('Parent ID do not match the current gene' + '\t'.join(line))
                    if strand == '+':
                        if int(line[4]) == end: line[4] = str(int(line[4]) + S)
                    if strand == '-':
                        if int(line[3]) == end: line[3] = str(max(0, int(line[3]) - S))
                fout.write('\t'.join(line) + '\n')








