import argparse

if __name__=='__main__':
    """
    Read output of bam-readcount and select all positions that have more than 3 non-reference bases mapped.
    
    The input file was generated using:
    bam-readcount -b 30 -q 20 -w 0 -f $refcur *.rmdup.sorted.bam | awk '{if($4>3) {n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":"); n6=split($11,f, ":"); n7=split($12,g, ":"); n8=split($13,h, ":");  print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2], f[1], f[2], g[1], g[2], h[1], h[2]}}' > bam_read_counts_b30_q20.txt
    
    The input file format is:
    Contig_ID   POS REF_Allele  READ_count  A_count C_count G_Count T_count N_Count  [INDEL_counts] [INDEL_counts]  [INDEL_counts]  
    """
    parser = argparse.ArgumentParser("Select positions which have more than three reads with non-reference alleles")
    parser.add_argument('f', help='path to bam_readcount output', type=argparse.FileType('r'))
    parser.add_argument('-n', dest='n', help='minimum number of non-reference reads', type=int, default=3)
    parser.add_argument('-p', dest='p', help='prefix for output file', type=str, default='filtered_low_ref_al_')
    parser.add_argument('-o', dest='o', help='Output file name', type=argparse.FileType('w'))

    args = parser.parse_args()
    from collections import deque

    with open(args.f.name) as fin:
        if args.o is None:
            if '/'.join(args.f.name.split("/")[:-1]) == '':
                foutname = './' + "{}".format(args.p) + args.f.name.split("/")[-1]
            else:
                foutname = '/'.join(args.f.name.split("/")[:-1]) + "/{}".format(args.p) + args.f.name.split("/")[-1]
        else:
            foutname = args.o.name
        with open(foutname, 'w') as fout:
            count = 0
            outstr = deque()
            base_dict = {
                "A": 4,
                "C": 5,
                "G": 6,
                "T": 7
            }
            for line in fin:
                l = line.strip().split(' ')

                ## Do not consider genomic positions which are N
                if l[2] == 'N':
                    continue

                if int(l[3]) - int(l[base_dict[l[2]]]) >= args.n:
                    outstr.append(line)
                    count += 1
                if count == 1000000:
                    fout.write(''.join(outstr))
                    outstr = deque()
                    count = 0
            fout.write(''.join(outstr))


