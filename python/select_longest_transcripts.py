def readfasta(f):
    from gzip import open as gzopen
    from collections import deque
    out = {}
    chrid = ''
    chrseq = deque()
    try:
        with gzopen(f, 'rb') as fin:
            for line in fin:
                break
                if b'>' in line:
                    if chrid != '':
                        out[chrid] = ''.join(chrseq)
                        chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                        chrseq = deque()
                    else:
                        chrid = line.strip().split(b'>')[1].split(b' ')[0].decode()
                else:
                    chrseq.append(line.strip().decode())
    except OSError:
        try:
            with open(f, 'r') as fin:
                for line in fin:
                    if '>' in line:
                        if chrid != '':
                            out[chrid] = ''.join(chrseq)
                            chrid = line.strip().split('>')[1].split(' ')[0]
                            chrseq = deque()
                        else:
                            chrid = line.strip().split('>')[1].split(' ')[0]
                    else:
                        chrseq.append(line.strip())
        except Exception as e:
            raise Exception(e)
    except Exception as e:
        raise Exception(e)

    if chrid != '':
        out[chrid] = ''.join(chrseq)
    # TODO: add check for the validation of input fasta files
    return out


def writefasta(fa, f):
    """
    :param fa: dictionary. Keys are chromosome ids. Values are sequence.
    :param f: Output file
    :return:
    """
    # TODO: Add capability to write fa.gzip files as well
    with open(f, 'w') as fo:
        for k, v in fa.items():
            fo.write('>'+k+'\n')
            for i in range(0, len(v), 60):
                fo.write(v[i:(i+60)]+'\n')

import argparse

if __name__ =='__main__':
    parser = argparse.ArgumentParser("Select longest transcript from genes with different isoforms", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', help='Fasta file containing isoform sequences. Different isoform ID should be named line Gene1.1, Gene1.2, Gene1.3 etc', type=argparse.FileType('r'))
    parser.add_argument('out', help='Output fasta file name', type=argparse.FileType('w'))

    args = parser.parse_args()

    fasta = readfasta(args.fasta.name)

    from collections import deque, defaultdict
    genelen = defaultdict(deque)
    for k,v in fasta.items():
        genelen['.'.join(k.split(".")[:-1])].append([k, len(v)])

    outfasta = dict()
    for k, v in genelen.items():
        if len(v) == 1: outfasta[v[0][0]] = fasta[v[0][0]]
        else:
            maxi = max(range(len(v)), key= lambda x: v[x][1])
            outfasta[v[maxi][0]] = fasta[v[maxi][0]]

    writefasta(outfasta, args.out.name)


