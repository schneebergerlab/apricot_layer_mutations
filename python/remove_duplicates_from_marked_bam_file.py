import argparse

"""
Take samtools markdup output (sorted by readname) and remove read-pairs for which both primary reads-alignments were found to be duplicated.
"""
if __name__ == '__main__':
    parser = argparse.ArgumentParser("Remove duplicated paired-end reads")
    parser.add_argument('bam', help='path to input bam file', type=argparse.FileType('r'))
    parser.add_argument('-p', dest='p', help='prefix for output fastq files', type=str, default='dedup')
    parser.add_argument('-c', dest='c', help='number of reads to print in a chunk', default=1000000, type=int)

    args = parser.parse_args()

    import pysam
    from collections import deque
    from gzip import open as gzopen

    b = pysam.AlignmentFile(args.bam.name, 'rb')
    # b = pysam.AlignmentFile('/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/tests/test_alignment_WT_1_using_bowtie2/AAACCTGAGGAATGGA.DUPmarked.deduped.bam', 'rb')
    # ls = [al for al in b if format(al.flag, '012b')[3] == '1']
    r1 = {}
    r2 = {}

    with open(args.p + "_R1.fastq", 'w') as f1:
        with open(args.p + "_R2.fastq", 'w') as f2:
            # with io.BufferedWriter(f1) as f1BF:
            count = 0
            r1_fq = deque()
            r2_fq = deque()
            for al in b:
                # print('here')
                # print(al.flag)
                if format(al.flag, '012b')[0] == '1' or format(al.flag, '012b')[3] == '1':
                    continue
                if al.is_read1:
                    r1[al.qname] = [al.seq if not al.is_reverse else al.get_forward_sequence(),
                                    al.qual if not al.is_reverse else al.qual[::-1],
                                    True if format(al.flag, '012b')[1] == '1' else False]
                elif al.is_read2:
                    r2[al.qname] = [al.seq if not al.is_reverse else al.get_forward_sequence(),
                                    al.qual if not al.is_reverse else al.qual[::-1],
                                    True if format(al.flag, '012b')[1] == '1' else False]
                else:
                    sys.exit("read is neither read1 or read2")
                if al.qname in r1.keys():
                    if al.qname in r2.keys():
                        if r1[al.qname][2] and r2[al.qname][2]:
                            r1.pop(al.qname)
                            r2.pop(al.qname)
                        else:
                            count += 1
                            r1_fq.extend(['@'+al.qname,
                                          r1[al.qname][0],
                                          '+',
                                          r1[al.qname][1]])
                            r2_fq.extend(['@'+al.qname,
                                          r2[al.qname][0],
                                          '+',
                                          r2[al.qname][1]])
                            r1.pop(al.qname)
                            r2.pop(al.qname)
                            if count % args.c == 0:
                                print(count, len(r1.keys()), len(r2.keys()))
                                f1.write('\n'.join(r1_fq) + '\n')
                                f2.write('\n'.join(r2_fq) + '\n')
                                r1_fq = deque()
                                r2_fq = deque()

            f1.write('\n'.join(r1_fq) + '\n')
            f2.write('\n'.join(r2_fq) + '\n')

# def test_write_fastq(seq):
#     with open('/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/tests/test.fa', 'w') as f:
#         f.write(seq)
#
# def test_write_fastq_gz(seq):
#     with gzopen('/srv/netscratch/dep_mercier/grp_schneeberger/projects/SynSearch/tests/test.fa.gz', 'wb') as f:
#         f.write(seq.encode())
#
# seq = '\n'.join(['>Chr1\n' + ''.join(np.random.choice(['A', 'T', 'G', 'C'], 150)) for i in range(1000)])
# %timeit test_write_fastq(seq)
# %timeit test_write_fastq_gz(seq)