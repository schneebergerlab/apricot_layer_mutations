import argparse

"""
Take samtools markdup output (sorted by readname) and remove read-pairs for which both primary reads-alignments were found to be duplicated.
"""
if __name__ == '__main__':
    parser = argparse.ArgumentParser("Remove duplicated paired-end reads")
    parser.add_argument('bam', help='path to input bam file', type=argparse.FileType('r'))
    parser.add_argument('-p', dest='p', help='prefix for output fastq files', type=str, default='dedup')

    args = parser.parse_args()

    import pysam
    from collections import deque
    from gzip import open as gzopen

    b = pysam.AlignmentFile(args.bam.name, 'rb')
    # ls = [al for al in b if format(al.flag, '012b')[3] == '1']
    r1 = {}
    r2 = {}

    with gzopen(args.p + "_R1.fastq.gz", 'wb') as f1:
        with gzopen(args.p + "_R2.fastq.gz", 'wb') as f2:
            count = 0
            r1_fq = deque()
            r2_fq = deque()
            for al in b:
                # print('here')
                if format(al.flag, '012b')[0] == '1' or format(al.flag, '012b')[3] == '1':
                    continue
                if al.is_read1:
                    r1[al.qname] = [al.seq, al.qual, True if format(al.flag, '012b')[1] == '1' else False]
                elif al.is_read2:
                    r2[al.qname] = [al.seq, al.qual, True if format(al.flag, '012b')[1] == '1' else False]
                else:
                    sys.exit("read is neither read1 or read2")
                if al.qname in r1.keys() and al.qname in r2.keys():
                    if r1[al.qname][2] and r2[al.qname][2]:
                        r1.pop(al.qname)
                        r2.pop(al.qname)
                    else:
                        count += 1
                        r1_fq.extend([b'@'+al.qname.encode(),
                                      r1[al.qname][0].encode(),
                                      b'+',
                                      r1[al.qname][1].encode()])
                        r2_fq.extend([b'@'+al.qname.encode(),
                                      r2[al.qname][0].encode(),
                                      b'+',
                                      r2[al.qname][1].encode()])
                if count == 100000:
                    count = 0
                    f1.write(b'\n'.join(r1_fq) + b'\n')
                    f2.write(b'\n'.join(r2_fq) + b'\n')
                    r1_fq = deque()
                    r2_fq = deque()

            f1.write(b'\n'.join(r1_fq) + b'\n')
            f2.write(b'\n'.join(r2_fq) + b'\n')
