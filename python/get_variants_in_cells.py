#!/usr/bin/env python
"""
Cluster the mutations (SNPs and indels) based on the barcode they are present in.
"""



import argparse

BASE_DICT = {
    "A": 4,
    "C": 5,
    "G": 6,
    "T": 7
}

def getrc(bcdir, smuts, mutfin, ref):
    from subprocess import Popen, PIPE
    from collections import deque
    bc = os.path.basename(bcdir)
    print(bc)
    p = Popen('bam-readcount -w 0 -b 0 -q 0 -l {mutfin} -f {ref} {bcdir}/{bc}.DUPmarked.deduped.bam'.format(mutfin=mutfin, ref=ref, bcdir=bcdir, bc=bc).split(), stdout=PIPE, stderr=PIPE)
    o = p.communicate()

    mutout = deque()
    for loc in o[0].decode().split('\n')[:-1]:
        loc = loc.split('\t')
        v = smuts[(loc[0], loc[1])][1]
        for j in loc[5:]:
            j = j.split(':')
            if j[0] != v: continue
            if int(j[1]) > 0: mutout.append(((loc[0], loc[1])))
    return {bc:mutout}
#end

if __name__ == '__main__':
    parser = argparse.ArgumentParser("Get allele depth of mutations in read-count data",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('mutations', help='List of mutations', type=argparse.FileType('r')) # BED-like file. For read position only the third column is used.
    parser.add_argument('ref', help='Reference genome', type=argparse.FileType('r'))
    parser.add_argument('-rc', help='bam-readcount output for the cell', type=argparse.FileType('r'))
    parser.add_argument('-frc', help='File containing paths to bam-readcount files. One per line', type=argparse.FileType('r'))
    # parser.add_argument('-n', dest='n', help='minimum alternate read count to select position', type=int, default=1)
    parser.add_argument('-o', dest='o', help='prefix for the output file', default='mut_rc', type=str)
    parser.add_argument('-p', dest='p', help='Number of CPU processors to use', default=1, type=int)

    args = parser.parse_args()

    import sys
    import os
    from collections import OrderedDict, defaultdict, deque, Counter
    from time import time
    from multiprocessing import Pool
    from functools import partial
    from tqdm import tqdm
    import pickle

    F1      = args.mutations.name
    # N       = args.n
    PRES    = args.o
    CWD = os.getcwd() + os.sep
    BCLIST = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/{}/barcodes_list'
    REF =

    if args.rc is not None and args.frc is not None:
        raise ValueError('-rc and -frc both cannot be provided. Using only -frc')
    if args.rc is None and args.frc is None:
        raise ValueError('Need bam-readcount path. Both -rc and -frc are not defined. Exiting.')
        sys.exit()
    if args.rc is not None: F2      = args.rc.name
    else : F2      = args.frc.name

    mlist = OrderedDict()
    with open(F1, 'r') as fin:
        # Data read for each mutation:
        # (chr, pos): (ref, alt, altcnt, af, sample)
        for line in fin:
            line = line.strip().split()
            mlist[(line[0], line[1])] = (line[2], line[3], int(line[4]), round(float(line[5]), 2), line[6])

    SAMPLES = set([v[4] for v in mlist.values()])
    for sample in SAMPLES:
        smuts = {k: v for k, v in mlist.items() if v[4] == sample}
        with open(BCLIST.format(sample), 'r') as fin:
            bcs = [line.strip() for line in fin]
        os.chdir(CWD)
        t = int(time())
        with open(str(t)+'TMP.bed', 'w') as fout:
            for k, v in smuts.items():
                fout.write("\t".join([k[0], k[1], k[1]]) + '\n')
        # with Pool(processes=args.p) as pool:
        with Pool(processes=3) as pool:
            # mutbc = tqdm(pool.map(partial(getrc, smuts=smuts, mutfin=str(t)+'TMP.bed', ref=args.ref.name), bcs))
            mutbc = pool.map(partial(getrc, smuts=smuts, mutfin=str(t)+'TMP.bed', ref=ref), tqdm(bcs))
        break
    # with open(CWD + 'mutbc.pickle', 'wb') as fout:
    #     pickle.dump(mutbc, fout)
    with open(CWD + 'mutbc.pickle', 'rb') as fout:
        mutbc = pickle.load(fout)


def getclusters(mutbc):
    import igraph as ig
    import random
    from functools import reduce
    from matplotlib import cm
    # Attempt one: Cluster the barcodes using mutations as features. This gives good looking clusters, however, mutations are not unique for the clusters. This could be a result of noise, but can also be a result of movement of cells across layers.
    muts = {k: v for i in mutbc for k,v in i.items()}
    b = sorted([k for k, v in muts.items() if len(v) > 0])
    G = ig.Graph()
    G.add_vertices(b)
    for i in range(len(b)):
        b1 = b[i]
        if len(muts[b1]) == 0: continue
        for j in range(i+1, len(b)):
            b2 = b[j]
            if len(muts[b2]) == 0: continue
            w = 0
            for k in muts[b1]:
                if k in muts[b2]: w+=1
            if w == 0: continue
            G.add_edge(b1, b2, weight=w)
    G.vs.select(_degree=0).delete()
    coords = G.layout('lgl')
    random.seed(1)
    G1 = G.community_leiden('modularity', weights='weight')
    cols = list(cm.Accent(range(len(set(G1.membership)))))
    ig.plot(G, vertex_color=[list(cols[i]) for i in G1.membership], edge_width=G.es['weight'], layout=coords).show()

    membership = G1.membership
    bcclust = defaultdict(deque)
    mutclust = defaultdict(deque)
    for i in range(len(b)):
        bcclust[membership[i]].append(b[i])
        mutclust[membership[i]].extend(muts[b[i]])
    for k, v in mutclust.items():
        print(k, Counter(v))

    allmut = reduce(lambda x,y: {*x, *y}, mutclust.values())

    from matplotlib import pyplot as plt
    y_off = [0]*len(allmut)
    for i in range(len(set(membership))):
        y = [mutclust[i].count(m) for m in allmut]
        plt.bar(range(len(allmut)), y, label="Number of BC: {}".format(len(bcclust[i])),bottom=y_off)
        y_off = [y[j] + y_off[j] for j in range(len(allmut))]
    plt.legend()
    mutanno = deque()
    for m in allmut:
        cs = [mutclust[i].count(m) for i in range(len(mutclust))]
        mutanno.append(cs.index(max(cs)))
    Counter(mutanno)
    # Counter({1: 7, 3: 3, 0: 3, 2: 2, 4: 1})


    # Attempt two: Cluster mutations using barcodes as features. This does not result in very good clusters as: 1) the number of mutations per sample is a quite low, 2) many mutations are supported by <10 barcodes
    # m = list(set([p for k, v in muts.items() for p in v]))
    muts2 = defaultdict(deque)
    for k, v in muts.items():
        for p in v: muts2['_'.join(p)].append(k)
    m = sorted(muts2.keys())
    G2 = ig.Graph()
    G2.add_vertices(m)
    for i in range(len(m)):
        m1 = m[i]
        for j in range(i+1, len(m)):
            m2 = m[j]
            w = 0
            for k in muts2[m1]:
                if k in muts2[m2]: w+=1
            if w == 0: continue
            G2.add_edge(m1, m2, weight=w)
    G2.vs.select(_degree=0).delete()
    coords = G2.layout('lgl')
    random.seed(1)
    G3 = G2.community_leiden('modularity', weights='weight')
    cols = list(cm.Accent(range(len(set(G3.membership)))))
    ig.plot(G2, vertex_color=[list(cols[i]) for i in G3.membership], edge_width=G2.es['weight'], layout=coords).show()







    mutclust[k] = set(v)
    allmut = reduce(lambda x,y: [*x, *y], mutclust.values())
    dups = [k for k, v in Counter(allmut). if v > 1]
    len(set(allmut))
    len(dups)

















    keys = list(muts.keys())
    if args.rc is not None:
        ref_rc, alt_rc = get_rc(F2, muts)
        with open(PRES + ".txt", 'w') as fout:
            fout.write('\t'.join(["chr", "start", "end", "ref", "alt", "depth", "freq"]) + '\t' + '\t'.join([F2+"_ref", F2+"_alt"]) + '\n')
            for i in range(len(keys)):
                k = keys[i].split('_')
                v = muts[keys[i]]
                fout.write('\t'.join([k[0], k[1], k[1], v[0], v[1], v[2], v[3], ref_rc[i], alt_rc[i]]) + '\n')
    else:
        t = ['ref', 'alt']
        with open(F2, 'r') as fin:
            rcs = OrderedDict()
            for line in fin:
                line = line.strip().split()
                rcs[line[1]] = get_rc(line[0], muts)
                print(line[1])
            samples = list(rcs.keys())

            with open(PRES + ".txt", 'w') as fout:
                fout.write('\t'.join(["chr", "start", "end", "ref", "alt", "depth", "freq"]) + '\t' + '\t'.join([s+';'+t[j] for s in samples for j in range(2)]) + '\n')
                for i in range(len(keys)):
                    k = keys[i].split('_')
                    v = muts[keys[i]]
                    fout.write('\t'.join([k[0], k[1], k[1], v[0], v[1], v[2], v[3]] + [rcs[s][j][i] for s in samples for j in range(2)]) + '\n')



    run("bam-readcount -w 0 -b 0 -q 0 -l ../indel_loci_for_readcount.txt -f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta {}.sorted.bt2.bam ".format(sample).split(), stdout=f)
    f.close()
    with open("TMP.txt", 'r') as fin:
        with open('all_good_readcount_indels.b0.q0.txt', 'w') as fout:
            for line in fin:
                line = line.strip().split()
                fout.write("\t".join(line[:4] + [i.split(':')[1] for i in line[5:10]] + [j for i in line[10:] for j in i.split(':')[:2] ]) + "\n")







vim

