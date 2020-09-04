"""
Take input paired input fasta files and filter out reads smaller than length L. Both pairs are filtered even if only one of them is smaller.
"""
import os
import sys
import argparse

parser = argparse.ArgumentParser('Filter small reads from the fastq files')
parser.add_argument('r1', help='R1 fasta file', type=argparse.FileType())
parser.add_argument('r2', help='R2 fasta file', type=argparse.FileType())
parser.add_argument('l', help='minimum read length', type=int)
parser.add_argument('-p', dest='p', help='prefix of output file', default=None)

args = parser.parse_args()

if args.p is None:
    args.p = 'L'+str(args.l)+'_'

from gzip import open as gzopen
from collections import deque

rrid1 = deque() ## remove read id

## Read the first read file and get read id of all small reads
c = 0
try:
    with gzopen(args.r1.name, 'rb') as fin:
        read = deque()
        for line in fin:
            c += 1
            if c % 1000000 == 0:
                print(c)
            read.append(line.strip())
            if len(read) == 4:
                rid, rs, s, rq = [i for i in read]                
                if len(rs) < args.l: rrid1.append(rid.split()[0])
                read.clear()
except OSError:
    with gzopen(args.r1.name, 'r') as fin:
        read = deque()
        for line in fin:
            read.append(line.strip())
            if len(read) == 4:
                rid, rs, s, rq = [i for i in read]
                if len(rs) < args.l: rrid1.append(rid.split()[0])
                read.clear()
except FileNotFoundError as e:
    raise e
'''
c = 0
try:
    with gzopen(args.r2.name, 'rb') as fin:
        read = deque()
        for line in fin:
            c += 1
            if c % 1000000 == 0:
                print(c)
            read.append(line.strip())
            if len(read) == 4:
                rid, rs, s, rq = [i for i in read]                
                if len(rs) < args.l: rrid.append(rid.split()[0])
                read.clear()
except OSError:
    with gzopen(args.r2.name, 'r') as fin:
        read = deque()
        for line in fin:
            read.append(line.strip())
            if len(read) == 4:
                rid, rs, s, rq = [i for i in read]
                if len(rs) < args.l: rrid.append(rid.split()[0])
                read.clear()
except FileNotFoundError as e:
    raise e

print(len(rrid))
sys.exit()
'''
rrid1 = set(list(rrid1))
## Read the second read file and get read id of all small reads and filter out matching reads from the file
rrid2 = deque()
try:
    with gzopen(args.r2.name, 'rb') as fin:
    # with gzopen('/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/mut/4747_C_merged_S9_L006_R2_001.fastq.gz', 'rb') as fin:
        with gzopen(args.p + os.path.basename(args.r2.name), 'w') as fout:
        # with gzopen('/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/mut/TEST_4747_C_merged_S9_L006_R2_001.fastq.gz', 'wb') as fout:
            read = deque()
            count = 0
            outstr = deque()
            c = 0
            for line in fin:
                c+=1
                if c % 1000000 == 0:
                    print(c)
                read.append(line.strip())
                if len(read) == 4:
                    rid, rs, s, rq = [i for i in read]
                    # if len(rs) >= args.l:
                    if len(rs) >= args.l:
                        if rid.split()[0] not in rrid1:
                            count += 1
                            outstr.append(rid+b'\n'+rs+b'\n'+s+b'\n'+rq+b'\n')
                    else: rrid2.append(rid.split()[0])
                    read.clear()
                    if count == 1000000:
                        count = 0
                        fout.write(b"".join(outstr))
                        outstr.clear()
            fout.write(b"".join(outstr))
except OSError:
    with open(args.r2.name, 'r') as fin:
        with open(args.p + os.path.basename(args.r2.name), 'w') as fout:
            read = deque()
            count = 0
            outstr = deque()
            for line in fin:
                read.append(line.strip())
                if len(read) == 4:
                    rid, rs, s, rq = [i for i in read]
                    # if len(rs) >= args.l:
                    if len(rs) >= args.l:
                        if rid.split()[0] not in rrid1:
                            count += 1
                            outstr.append(rid+'\n'+rs+'\n'+s+'\n'+rq+'\n')
                    else: rrid2.append(rid.split()[0])
                    read.clear()
                    if count == 1000000:
                        count = 0
                        fout.write("".join(outstr))
                        outstr.clear()
            fout.write("".join(outstr))
except FileNotFoundError as e:
    raise e

rrid1 = set(list(rrid2) + list(rrid1))

## Read the first read file and remove reads corresponding to small reads from both files
try:
    with gzopen(args.r1.name, 'rb') as fin:
    # with gzopen('/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/mut/4747_C_merged_S9_L006_R2_001.fastq.gz', 'rb') as fin:
        with gzopen(args.p + os.path.basename(args.r1.name), 'w') as fout:
        # with gzopen('/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/mut/TEST_4747_C_merged_S9_L006_R2_001.fastq.gz', 'wb') as fout:
            read = deque()
            count = 0
            outstr = deque()
            c = 0
            for line in fin:
                c += 1
                if c % 1000000 == 0:
                    print(c)
                read.append(line.strip())
                if len(read) == 4:
                    rid, rs, s, rq = [i for i in read]
                    # if len(rs) >= args.l:
                    # if len(rs) >= 50:
                    if rid.split()[0] not in rrid1:
                        count += 1
                        outstr.append(rid+b'\n'+rs+b'\n'+s+b'\n'+rq+b'\n')
                    # else: rrid2.append(rid.split()[0])
                    read.clear()
                    if count == 1000000:
                        count = 0
                        fout.write(b"".join(outstr))
                        outstr.clear()
            fout.write(b"".join(outstr))
except OSError:
    with open(args.r1.name, 'r') as fin:
        with open(args.p + os.path.basename(args.r1.name), 'w') as fout:
            read = deque()
            count = 0
            outstr = deque()
            for line in fin:
                read.append(line.strip())
                if len(read) == 4:
                    rid, rs, s, rq = [i for i in read]
                    # if len(rs) >= args.l:
                    # if len(rs) >= 50:
                    if rid.split()[0] not in rrid1:
                        count += 1
                        outstr.append(rid+'\n'+rs+'\n'+s+'\n'+rq+'\n')
                    # else: rrid2.append(rid.split()[0])
                    read.clear()
                    if count == 1000000:
                        count = 0
                        fout.write("".join(outstr))
                        outstr.clear()
            fout.write("".join(outstr))
except FileNotFoundError as e:
    raise e

