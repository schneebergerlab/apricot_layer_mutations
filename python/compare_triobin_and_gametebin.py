import argparse
from collections import deque
import gzip

QS= "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"

TRIO_CUR = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/test/canu_trio/canu_trio/haplotype/haplotype-cur.fasta.gz'
TRIO_ORA = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/test/canu_trio/canu_trio/haplotype/haplotype-ora.fasta.gz'
TRIO_UNK = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/test/canu_trio/canu_trio/haplotype/haplotype-unknown.fasta.gz'


def getreadids(f, gz=False):
    rids = deque()
    if gz:
        with gzip.open(f, 'rb') as fin:
            for line in fin:
                line = line.strip()
                if line[0] == b'>'[0]:
                    line = line.decode()
                    try:
                        end = line.index(' ')
                    except ValueError as e:
                        end = len(line)
                    rids.append(line[1:end])
    else:
        with open(f, 'r') as fin:
            for line in fin:
                line = line.strip()
                if line[0] == '>':
                    try:
                        end = line.index(' ')
                    except ValueError as e:
                        end = len(line)
                    rids.append(line[1:end])
    return set(rids)

if __name__ == '__main__':
    """
    Compare the binning of reads from trio-binning and gamete-binning.
    
    """
    parser = argparse.ArgumentParser("Select positions which are supported by read-level features",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('reg', help='List of regions to check', type=argparse.FileType('r'))
    parser.add_argument('bam', help='bam file', type=argparse.FileType('r'))
    parser.add_argument('ref', help='path to reference file', type=argparse.FileType('r'))

    parser.add_argument('-n', dest='n', help='minimum alternate read count to select position', type=int, default=1)
    parser.add_argument('-o', dest='o', help='prefix for the output file', default='cell_var', type=str)

    args = parser.parse_args()

    trio_cur_rids = getreadids(TRIO_CUR, True)
    trio_ora_rids = getreadids(TRIO_ORA, True)
    trio_unk_rids = getreadids(TRIO_UNK, True)

    rid_dist = {}
    indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/gamete_binning/apricot_hifi_snp_marker_separated_pbreads/'
    for i in range(1, 9):
        for c in ['PPP', 'MMM']:
            sample = str(i) + '.txt_' + c + '_pbreads.fa'
            print(sample)
            rids = getreadids(indir + sample)
            cur = rids.intersection(trio_cur_rids)
            ora = rids.intersection(trio_ora_rids)
            unk = rids.intersection(trio_unk_rids)
            rid_dist[sample] = (len(cur), len(ora), len(unk))

    for sample in ['mapped_no_lg_pb_reads.fa', 'unmapped_starcigar_pb_reads.fa']:
        rids = getreadids(indir + sample)
        cur = rids.intersection(trio_cur_rids)
        ora = rids.intersection(trio_ora_rids)
        unk = rids.intersection(trio_unk_rids)
        rid_dist[sample] = (len(cur), len(ora), len(unk))

    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=[8, 10])
    count = 0
    for i in range(1, 9):
        for c in ['PPP', 'MMM']:
            count += 1
            sample = str(i) + '.txt_' + c + '_pbreads.fa'
            ax = fig.add_subplot(5,  4, count)
            t = sum(rid_dist[sample])
            ax.bar(['cur', 'ora', 'unk'],
                   [v/t for v in rid_dist[sample]])
            ax.set_xlabel(sample)

    for sample in ['mapped_no_lg_pb_reads.fa', 'unmapped_starcigar_pb_reads.fa']:
        count += 1
        ax = fig.add_subplot(5,  4, count)
        t = sum(rid_dist[sample])
        ax.bar(['cur', 'ora', 'unk'],
               [v/t for v in rid_dist[sample]])
        ax.set_xlabel(sample)
        count += 1

    plt.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.98, wspace=0.5, hspace=0.5)
    plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/gametebinning_vs_triobinning.pdf')

    cur_lg = deque()
    ora_lg = deque()
    trio_cur = [0, 0, 0, 0]
    trio_ora = [0, 0, 0, 0]
    for i in range(1, 9):
        for c in ['PPP', 'MMM']:
            sample = str(i) + '.txt_' + c + '_pbreads.fa'
            if rid_dist[sample][0] > rid_dist[sample][1]:
                cur_lg.append(sample)
                trio_cur[0] += rid_dist[sample][0]
                trio_ora[1] += rid_dist[sample][1]
            else:
                ora_lg.append(sample)
                trio_cur[1] += rid_dist[sample][0]
                trio_ora[0] += rid_dist[sample][1]
    trio_cur[2] = rid_dist['mapped_no_lg_pb_reads.fa'][0]
    trio_cur[3] = rid_dist['unmapped_starcigar_pb_reads.fa'][0]

    trio_ora[2] = rid_dist['mapped_no_lg_pb_reads.fa'][1]
    trio_ora[3] = rid_dist['unmapped_starcigar_pb_reads.fa'][1]



