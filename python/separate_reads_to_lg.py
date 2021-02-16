import pysam as sam
from collections import deque
import os
from datetime import datetime

curdir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/contig_grouping/cur/'
oradir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/contig_grouping/ora/'
curlgs = {}
oralgs = {}
for i in range(1, 9):
    with open("{}final_manual_upd_map_group{}_contigs.fasta".format(curdir, i), 'r') as fin:
        contigs = deque()
        for line in fin:
            if line[0] != '>': continue
            line = line.strip()
            contigs.append(line[1:])
        for contig in contigs:
            curlgs[contig] = i

    with open("{}final_manual_upd_map_group{}_contigs.fasta".format(oradir, i), 'r') as fin:
        contigs = deque()
        for line in fin:
            if line[0] != '>': continue
            line = line.strip()
            contigs.append(line[1:])
        for contig in contigs:
            oralgs[contig] = i
contigs_list = set(list(curlgs.keys()) + list(oralgs.keys()))

os.makedirs("curlgs", exist_ok=True)
os.makedirs("oralgs", exist_ok=True)
curfoutR1 = []
curfoutR2 = []
orafoutR1 = []
orafoutR2 = []
curfoutR1str = [deque() for i in range(8)]
curfoutR2str = [deque() for i in range(8)]
orafoutR1str = [deque() for i in range(8)]
orafoutR2str = [deque() for i in range(8)]
for i in range(1, 9):
    os.makedirs("curlgs/lg{}".format(i), exist_ok=True)
    curfoutR1.append(open('curlgs/lg{}/R1.fq'.format(i), 'w'))
    curfoutR2.append(open('curlgs/lg{}/R2.fq'.format(i), 'w'))
    os.makedirs("oralgs/lg{}".format(i), exist_ok=True)
    orafoutR1.append(open('oralgs/lg{}/R1.fq'.format(i), 'w'))
    orafoutR2.append(open('oralgs/lg{}/R2.fq'.format(i), 'w'))

bamfile = sam.AlignmentFile('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/lg_read_set/rp_diploid.hic.n_sorted.bam', 'rb')
# bamfile = sam.AlignmentFile('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/lg_read_set/TMP.sorted.bam', 'rb')
cnt = 0
qname = ''
R1seq = ''
R2seq = ''
R1qual = ''
R2qual = ''
readcurlgs = deque()
readoralgs = deque()
cnt2 = 0
for read in bamfile.fetch(until_eof=True):
    cnt+=1
    if cnt == 1000000:
        cnt2+=cnt
        print(cnt2, str(datetime.now()))
        for i in range(8):
            if len(curfoutR1str[i]) > 0:
                curfoutR1[i].write('\n'.join(curfoutR1str[i]) + '\n')
            if len(curfoutR2str[i]) > 0:
                curfoutR2[i].write('\n'.join(curfoutR2str[i]) + '\n')
            if len(orafoutR1str[i]) > 0:
                orafoutR1[i].write('\n'.join(orafoutR1str[i]) + '\n')
            if len(orafoutR2str[i]) > 0:
                orafoutR2[i].write('\n'.join(orafoutR2str[i]) + '\n')
            curfoutR1str[i] = deque()
            curfoutR2str[i] = deque()
            orafoutR1str[i] = deque()
            orafoutR2str[i] = deque()
        cnt = 0

    if read.qname != qname:
        if qname == '':
            qname = read.qname
            flag = format(read.flag, '012b')
            if flag[0] != '1' and flag[2] != '1':
                if read.reference_name in contigs_list:
                    if 'cur' in read.reference_name:
                        readcurlgs.append(curlgs[read.reference_name])
                    else:
                        readoralgs.append(oralgs[read.reference_name])
            if read.seq is not None:
                if flag[5] == '1':
                    R1seq = read.get_forward_sequence()
                    R1qual = read.qual[::-1] if read.is_reverse else read.qual
                else:
                    R2seq = read.get_forward_sequence()
                    R2qual = read.qual[::-1] if read.is_reverse else read.qual
        else:
            for lgs in set(readcurlgs):
                curfoutR1str[lgs-1].extend(['@' + qname, R1seq, '+', R1qual])
                curfoutR2str[lgs-1].extend(['@' + qname, R2seq, '+', R2qual])
            for lgs in set(readoralgs):
                orafoutR1str[lgs-1].extend(['@' + qname, R1seq, '+', R1qual])
                orafoutR2str[lgs-1].extend(['@' + qname, R2seq, '+', R2qual])
            readcurlgs = deque()
            readoralgs = deque()
            qname = read.qname
            flag = format(read.flag, '012b')
            if flag[0] != '1' and flag[2] != '1':
                if read.reference_name in contigs_list:
                    if 'cur' in read.reference_name:
                        readcurlgs.append(curlgs[read.reference_name])
                    else:
                        readoralgs.append(oralgs[read.reference_name])
            if read.seq is not None:
                if flag[5] == '1':
                    R1seq = read.get_forward_sequence()
                    R1qual = read.qual[::-1] if read.is_reverse else read.qual
                else:
                    R2seq = read.get_forward_sequence()
                    R2qual = read.qual[::-1] if read.is_reverse else read.qual
    else:
        flag = format(read.flag, '012b')
        if flag[0] != '1' and flag[2] != '1':
            if read.reference_name in contigs_list:
                if 'cur' in read.reference_name:
                    readcurlgs.append(curlgs[read.reference_name])
                else:
                    readoralgs.append(oralgs[read.reference_name])
        if read.seq is not None:
            if flag[5] == '1':
                R1seq = read.get_forward_sequence()
                R1qual = read.qual[::-1] if read.is_reverse else read.qual
            else:
                R2seq = read.get_forward_sequence()
                R2qual = read.qual[::-1] if read.is_reverse else read.qual

for lgs in set(readcurlgs):
    curfoutR1str[lgs-1].extend(['@' + qname, R1seq, '+', R1qual])
    curfoutR2str[lgs-1].extend(['@' + qname, R2seq, '+', R2qual])
for lgs in set(readoralgs):
    orafoutR1str[lgs-1].extend(['@' + qname, R1seq, '+', R1qual])
    orafoutR2str[lgs-1].extend(['@' + qname, R2seq, '+', R2qual])
if len(curfoutR1str[i]) > 0:
    curfoutR1[i].write('\n'.join(curfoutR1str[i]) + '\n')
if len(curfoutR2str[i]) > 0:
    curfoutR2[i].write('\n'.join(curfoutR2str[i]) + '\n')
if len(orafoutR1str[i]) > 0:
    orafoutR1[i].write('\n'.join(orafoutR1str[i]) + '\n')
if len(orafoutR2str[i]) > 0:
    orafoutR2[i].write('\n'.join(orafoutR2str[i]) + '\n')

[f.close() for f in curfoutR1]
[f.close() for f in curfoutR2]
[f.close() for f in orafoutR1]
[f.close() for f in orafoutR2]