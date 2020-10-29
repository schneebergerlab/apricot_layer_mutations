import argparse

if __name__=='__main__':
    parser=argparse.ArgumentParser("Trim 10x barcodes from the 5' end of the R1 reads")
    parser.add_argument('bc', help='length of the barcode', type=int)
    parser.add_argument('f1', help='path to R1_fastq.gz file', type=argparse.FileType('r'))
    parser.add_argument('f2', help='path to R2_fastq.gz file', type=argparse.FileType('r'))
    parser.add_argument('-p', dest='p', help='prefix for output file', type=str, default='bc_trimmed_')

    args = parser.parse_args()

    from gzip import open as gzopen
    from collections import deque

    with gzopen(args.f1.name, 'rb') as fin:
        if '/'.join(args.f1.name.split("/")[:-1]) == '':
            foutname = './' + "{}".format(args.p) + args.f1.name.split("/")[-1]
        else:
            foutname = '/'.join(args.f1.name.split("/")[:-1]) + "/{}".format(args.p) + args.f1.name.split("/")[-1]
        with gzopen(foutname, 'wb') as fout:
            count=0
            outstr = deque()
            for line in fin:
                count = (count+1) % 4
                if count == 1:
                    outstr.append(line)
                elif count == 2:
                    outstr.append(line[args.bc:])
                elif count == 3:
                    outstr.append(line)
                elif count == 0:
                    outstr.append(line[args.bc:])
                if len(outstr) > 1000000:
                    fout.write(b''.join(outstr))
                    outstr = deque()
            fout.write(b''.join(outstr))

    import shutil
    if '/'.join(args.f2.name.split("/")[:-1]) == '':
        foutname = './' + "{}".format(args.p) + args.f2.name.split("/")[-1]
    else:
        foutname = '/'.join(args.f2.name.split("/")[:-1]) + "/{}".format(args.p) + args.f2.name.split("/")[-1]

    shutil.copyfile(args.f2.name, foutname)