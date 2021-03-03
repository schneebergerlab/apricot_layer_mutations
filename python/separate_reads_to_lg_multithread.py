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

def get_contig_reads(bam, contig):
    import pysam
    bamfile = sam.AlignmentFile('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/lg_read_set/rp_diploid.hic.sorted.bam', 'rb')

bam = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/lg_read_set/rp_diploid.hic.n_sorted.bam'

