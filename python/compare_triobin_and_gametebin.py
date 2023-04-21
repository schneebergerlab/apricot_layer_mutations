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



#### MERYL and KMC based KMER validation for gamete-binning and trio-binning kmers

# Meryl results for gamete-binning
INDIR='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/gamete_binning/apricot_hifi_snp_marker_separated_pbreads/'

meryl_gb = {}
for i in range(1, 9):
    for s in ['PPP', 'MMM']:
        cursum = 0
        orasum = 0
        with open(INDIR + str(i) + '_' + s + '_hapmers/currot.inherited.hist') as f:
            for line in f:
                line = line.strip().split()
                if int(line[0]) > 5:
                    cursum += int(line[1])
        with open(INDIR + str(i) + '_' + s + '_hapmers/orangered.inherited.hist') as f:
            for line in f:
                line = line.strip().split()
                if int(line[0]) > 5:
                    orasum += int(line[1])
        meryl_gb[str(i) + '_' + s] = [cursum/(cursum+orasum), orasum/(cursum+orasum)]
        print(cursum, orasum)

# KMC results for gamete-binning
kmc_gb = {}
for i in range(1, 9):
    for s in ['PPP', 'MMM']:
        cursum = 0
        orasum = 0
        with open('{indir}/intersect_{i}_{s}/{i}_{s}_currot_21mers.histo'.format(indir=INDIR, i=i, s=s)) as f:
            for line in f:
                line = line.strip().split()
                if int(line[0]) > 5:
                    cursum += int(line[1])
        with open('{indir}/intersect_{i}_{s}/{i}_{s}_orangered_21mers.histo'.format(indir=INDIR, i=i, s=s)) as f:
            for line in f:
                line = line.strip().split()
                if int(line[0]) > 5:
                    orasum += int(line[1])

        kmc_gb[str(i) + '_' + s] = [cursum/(cursum+orasum), orasum/(cursum+orasum)]
        print(cursum, orasum)

# Meryl results for trio-binning
INDIR='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/trio_binning/'

meryl_tb = {}
for s in ['cur', 'ora']:
    cursum = 0
    orasum = 0
    with open('{indir}/{s}_hapmers/currot.inherited.hist'.format(indir=INDIR, s=s)) as f:
        for line in f:
            line = line.strip().split()
            if int(line[0]) > 5:
                cursum += int(line[1])

    with open('{indir}/{s}_hapmers/orangered.inherited.hist'.format(indir=INDIR, s=s)) as f:
        for line in f:
            line = line.strip().split()
            if int(line[0]) > 5:
                orasum += int(line[1])
    meryl_tb[s] = [cursum/(cursum+orasum), orasum/(cursum+orasum)]

# KMC results for trio-binning
kmc_tb = {}
for s in ['cur', 'ora']:
    cursum = 0
    orasum = 0
    with open('{indir}/intersect_{s}/{s}_currot_21mers.histo'.format(indir=INDIR, s=s)) as f:
        for line in f:
            line = line.strip().split()
            if int(line[0]) > 5:
                cursum += int(line[1])

    with open('{indir}/intersect_{s}/{s}_orangered_21mers.histo'.format(indir=INDIR, s=s)) as f:
        for line in f:
            line = line.strip().split()
            if int(line[0]) > 5:
                orasum += int(line[1])
    kmc_tb[s] = [cursum/(cursum+orasum), orasum/(cursum+orasum)]

markers = {1:"v",
           2:"^",
           3:"<",
           4:">",
           5:"1",
           6:"2",
           7:"3",
           8:"4",
           }
colors = {"PPP" : 'red',
          "MMM" : 'black'
          }

from matplotlib import pyplot as plt
fig = plt.figure(figsize=[8, 10])
# fig = plt.figure()

ax = fig.add_subplot(2, 2, 1)
for i in range(1, 9):
    for s in ['PPP', 'MMM']:
        ax.scatter([meryl_gb[str(i)+"_"+s][0]],
                   [meryl_gb[str(i)+"_"+s][1]],
                   c=colors[s],
                   marker=markers[i],
                   label=str(i)+"_"+s)
ax.set_xlim([0,1])
ax.legend(ncol=2)
ax.set_xlabel("parent cur")
ax.set_ylabel("parent ora")
ax.set_title("Meryl validation of GB phased reads")

ax = fig.add_subplot(2, 2, 3)
for i in range(1, 9):
    for s in ['PPP', 'MMM']:
        ax.scatter([kmc_gb[str(i)+"_"+s][0]],
                   [kmc_gb[str(i)+"_"+s][1]],
                   c=colors[s],
                   marker=markers[i],
                   label=str(i)+"_"+s)
ax.set_xlim([0,1])
ax.legend(ncol=2)
ax.set_xlabel("parent cur")
ax.set_ylabel("parent ora")
ax.set_title("KMC validation of GB phased reads")

ax = fig.add_subplot(2, 2, 2)
ax.scatter([v[0] for v in meryl_tb.values()],
           [v[1] for v in meryl_tb.values()],
           c = "black")
ax.set_xlabel("parent cur")
ax.set_ylabel("parent ora")
ax.set_title("Meryl validation of TB phased reads")


ax = fig.add_subplot(2, 2, 4)
ax.scatter([v[0] for v in kmc_tb.values()],
           [v[1] for v in kmc_tb.values()],
           c="black")
ax.set_xlabel("parent cur")
ax.set_ylabel("parent ora")
ax.set_title("KMC validation of TB phased reads")

plt.subplots_adjust(left=0.075, bottom=0.05, right=0.98, top=0.95, wspace=0.2, hspace=0.2)
plt.savefig('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/meryl_kmc_validation_of_gb_and_tb_phasing.pdf')

indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/read_phasing/gamete_binning/apricot_hifi_snp_marker_separated_pbreads/'
for i in range(1, 9):
    for s in ['PPP', 'MMM']:
        with open("{indir}/intersect_{i}_{s}/{i}_{s}_currot_21mers.count".format(indir=indir, i=i, s=s)) as f:
            for line in f:
                cur_cnt = int(line.strip().split()[0])

        with open("{indir}/intersect_{i}_{s}/{i}_{s}_orangered_21mers.count".format(indir=indir, i=i, s=s)) as f:
            for line in f:
                ora_cnt = int(line.strip().split()[0])
        print(str(i)+s, cur_cnt, ora_cnt, (cur_cnt/ (cur_cnt + ora_cnt)), (ora_cnt/ (cur_cnt + ora_cnt)))
