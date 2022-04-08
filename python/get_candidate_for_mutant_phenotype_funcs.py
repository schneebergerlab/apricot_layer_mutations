import argparse

def get_pos(fin):
    from collections import defaultdict, deque
    with open(fin, 'r') as f:
        out = defaultdict(deque)
        count = 0
        for line in f:
            count += 1
            line = line.strip().split()
            out[line[0]].append(line[1])
    for k, v in out.items():
        out[k] = set(v)
    return out
#END

def write_pos_only_snps(fin, fout, pos, CAN_N, RC_MIN, RC_MAX, CAN_CNT=-1):
    '''
    Select all positions that one alternate substitution allele (SNPs).
    '''
    import pandas as pd
    from collections import deque
    from tqdm import tqdm
    base_dict = {4: 'A',
                 5: 'C',
                 6: 'G',
                 7: 'T'
                 }
    
    with open(fin, 'r') as f:
        df = deque()
        for line in tqdm(f):
            line = line.strip().split()
            if int(line[3]) < RC_MIN: continue
            if int(line[3]) > RC_MAX: continue
            hf = 0
            # for i in [4, 5, 6, 7, 10, 12, 14]:
            for i in [4, 5, 6, 7]:
                try:
                    if int(line[i]) >= CAN_N:
                        hf += 1
                except IndexError as e:
                    break
            if hf > 2:
                continue
            try:
                if line[1] in pos[line[0]]:
                    ## Get alternate allele readcount and alternate allele ratio
                    rc = 0
                    ar = 0
                    alt = ''
                    for i in [4, 5, 6, 7]:
                        try:
                            if int(line[i]) >= CAN_N:
                                if base_dict[i] != line[2]:
                                    rc = int(line[i])
                                    ar = int(line[i])/int(line[3])
                                    alt = base_dict[i]
                        except IndexError as e:
                            pass
                    if rc == 0 or ar == 0 or alt == '':
                        continue
                    df.append([line[0], line[1], line[1], line[2], alt, rc, ar])
            except KeyError as e:
                pass

    df = pd.DataFrame(df)
    df[1] = df[1].astype(int)
    df[2] = df[2].astype(int)
    df[5] = df[5].astype(int)
    df_cand = df.iloc[0:CAN_CNT] if CAN_CNT != -1 else df.copy()
    df_cand.to_csv(fout+'.regions', sep=' ', index=False, header=False)
    df_bed = pd.DataFrame({'chr'  : df_cand[0],
                           'start': df_cand[1]-1,
                           'end'  : df_cand[2],
                           'ref'  : df_cand[3],
                           'alt'  : df_cand[4],
                           'vaf'  : df_cand[5],
                           'var'  : df_cand[6]})
    # df_bed.to_csv(foutname+'.bed', sep='\t', index=False, header=False)
    df_bed.sort_values(['chr', 'start', 'end'], inplace=True)
    df_bed.to_csv(fout +'.sorted.bed', sep='\t', index=False, header=False)
#END

def sample_snps(fin, fout, CAN_N, RC_MIN, RC_MAX, CAN_CNT=-1):
    """
    fin = input file
    fout = output file name
    CAN_N = minimum alt read count
    RC_MIN = Minimum coverage at position
    RC_MAX = Maximum coverage at position
    CAN_CNT = Number of candidiates to output
    """
    pos = get_pos(fin)
    write_pos_only_snps(fin, fout, pos, CAN_N, RC_MIN, RC_MAX, CAN_CNT)
#END


def getlocs(f, n=20):
    '''
    Read bamrc or filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt file to select positions having more than n alt reads
    '''
    from collections import deque
    locs = deque()
    BASE_DICT = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
    with open(f, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            if int(line[3]) - int(line[BASE_DICT[line[2]]]) >= n:
                locs.append(f'{line[0]}_{line[1]}')
            elif len(line) > 9:
                RS = sum(map(int, line[4:8])) - int(line[BASE_DICT[line[2]]])
                RS += sum([int(line[i]) for i in range(10, len(line), 2)])
                if RS >= n:
                    locs.append(f'{line[0]}_{line[1]}')
    return set(locs)
#END


def filterbg(fg, bg, bgtype='any'):
    '''
    Select positions present in all foreground list of positions, and then
    filter out positions in bg.

    bgtype can be 'any' if a position if to be filtered out if it is present in any bg sample or 'all' if a position has to be filtered out only when it is present in all of the bg samples.
    '''
    fgmerge = fg[0].copy()
    for i in range(1, len(fg)):
        fgmerge = fgmerge.intersection(fg[i])
    bgmerge = bg[0]
    if bgtype == 'any':
        for i in range(1, len(bg)):
            bgmerge = bgmerge.union(bg[i])
    elif bgtype == 'all':
        for i in range(1, len(bg)):
            bgmerge = bgmerge.intersection(bg[i])
    else:
        print(f'Invalid bgtype value "{bgtype}". Only "any" and "all" are accepted.')
    fgmerge -= bgmerge
    return fgmerge
# END


def filterbgcnt(fg, bg, cnt=1):
    '''
    Select positions present in all foreground list of positions, and then
    filter out positions in bg.

    If the position is present in more than cnt samples, then it would be filtered out.
    '''
    from collections import Counter, deque
    fgmerge = fg[0].copy()
    for i in range(1, len(fg)):
        fgmerge = fgmerge.intersection(fg[i])
    bgs = list(bg[0])
    for i in range(1, len(bg)):
        bgs += list(bg[i])
    bgs = Counter(bgs)
    outfg = deque()
    for pos in fgmerge:
        try:
            c = bgs[pos]
        except KeyError:
            outfg.append(pos)
            continue
        if c > cnt:
            continue
        outfg.append(pos)
    return set(outfg)
# END

def filterclosepos(pos):
    from pandas import DataFrame as pd_df
    from pandas import concat
    p = [i.rsplit('_', 1) for i in pos]
    p = pd_df(p)
    p[1] = p[1].astype(int)
    pout = pd_df()
    for g in p.groupby(0):
        gc = g[1].copy()
        gc.sort_values([1], inplace=True)
        gc.reset_index(inplace=True, drop=True)
        g1 = [False] + list(gc[1].array[1:] - gc[1].array[:-1] < 100)
        g2 = list(gc[1].array[:-1] - gc[1].array[1:] > -100) + [False]
        gc['filter'] = [a or b for a, b in zip(g1, g2)]
        gc = gc.loc[gc['filter'] != True]
        pout = concat([pout, gc])
    return set(pout[0].astype(str) + '_' + pout[1].astype(str))
# END


def plot_selected_pos(pos, igvb, outdir, M=200, HEIGHT=600, emptydir=False):
    from subprocess import run
    import os
    indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
    # BAMS = ['WT_1/WT_1.sorted.bt2.bam',
    #         'wt7/wt7.deduped.bam',
    #         'wt18/wt18.deduped.bam',
    #         'WT_19/WT_19.sorted.bt2.bam',
    #         'mut4/mut4.deduped.bam',
    #         'MUT_11_1/MUT_11_1.sorted.bt2.bam',
    #         'mut11_2/mut11_2.deduped.bam',
    #         'MUT_15/MUT_15.sorted.bt2.bam']
    BAMS = ['wt7/wt7.deduped.bam',
            'wt18/wt18.deduped.bam',
            'mut4/mut4.deduped.bam',
            'mut11_2/mut11_2.deduped.bam']

    GENOME = "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta"
    IGV = '/srv/netscratch/dep_mercier/grp_schneeberger/software/igv/IGV_Linux_2.10.2/igv.sh'
    OUTDIR = outdir
    # Delete the contents of the output folder
    if emptydir:
        fins = os.listdir(outdir)
        for f in fins:
            os.remove(f'{os.getcwd()}/{outdir}/{f}')

    with open(igvb, 'w') as fout:
        fout.write("new"+"\n")
        fout.write("genome {}\n".format(GENOME))
        for bam in BAMS:
            fout.write("load {}\n".format(bam))
        fout.write("snapshotDirectory {}\n".format(OUTDIR))
        fout.write("maxPanelHeight {}\n".format(HEIGHT))
        fout.write("squish\n")
        for p in pos:
            c, m = p.rsplit('_', 1)
            s = int(m) - M
            e = int(m) + M
            fout.write("goto {}:{}-{}\n".format(c, s, e))
            fout.write("snapshot {}:{}-{}.png\n".format(c, s, e))
        fout.write("exit\n")
    with open(indir+"igv.log", 'a') as igvlog:
        run("{} -b {}".format(IGV, igvb).split(), stdout=igvlog, stderr=igvlog)
#END


def getlinefrombamrc(f, cand):
    from collections import deque
    from tqdm import tqdm
    out = deque()
    with open(f, 'r') as fin:
        for line in tqdm(fin):
            line = line.strip().split()
            if line[0] + "_" + line[1] not in cand:
                continue
            out.append(line)
    return out
#END

def generate_TEfiles_splitreader(tegff, famout):
    from collections import defaultdict, deque
    # tsdfin = '/srv/netscratch/dep_mercier/grp_schneeberger/software/splitreader/thaliana/superfamily_TSD.txt'
    # tegff = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/TE_annotation_edited.gff'
    # tsds = {}
    # with open(tsdfin, 'r') as fin:
    #     line = fin.readline()
    #     for line in fin:
    #         line = line.strip().split()
    #         tsds[line[0]] = int(line[1])
    gff = defaultdict(deque)
    with open(tegff, 'r') as fin:
        for i in range(3): fin.readline()
        for line in fin:
            line = line.strip().split()
            d = line[8].split(';')
            gff[d[2].split("=")[1] + '_' + d[3].split("=")[1]].append(d[1].split("=")[1])
    # with open("/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/TEfamily_superfamily.txt", 'w') as fout:
    with open(famout, 'w') as fout:
        for k,v in gff.items():
            for f in set(v):
                fout.write("\t".join([f, k]) + "\n")
#END


def tepid_te_del_igv_batch(bed, igvb, outdir):
    from collections import deque
    from subprocess import run
    indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/'
    BAMS = deque()
    M = 200
    for s in ('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15'):
        BAMS.append(indir + '{s}/{s}.bam'.format(s=s))
        # BAMS.append(indir + '{s}/{s}.discordants.sorted.bam'.format(s=s))
    GENOME = "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta"
    HEIGHT = 600
    IGV = '/srv/netscratch/dep_mercier/grp_schneeberger/software/igv/IGV_Linux_2.10.2/igv.sh'
    OUTDIR = outdir

    with open(bed, 'r') as fin, open(igvb, 'w') as fout:
        fout.write("new"+"\n")
        fout.write("genome {}\n".format(GENOME))
        for bam in BAMS:
            fout.write("load {}\n".format(bam))
        fout.write("snapshotDirectory {}\n".format(OUTDIR))
        fout.write("maxPanelHeight {}\n".format(HEIGHT))
        fout.write("squish\n")
        for line in fin:
            line = line.strip().split()
            c = line[0]
            s = int(line[1]) - M
            e = int(line[2]) + M
            fout.write("goto {}:{}-{}\n".format(c, s, e))
            fout.write("snapshot {}:{}-{}.png\n".format(c, s, e))
        fout.write("exit\n")
    with open(indir+"igv.log", 'a') as igvlog:
        run("{} -b {}".format(IGV, igvb).split(), stdout=igvlog, stderr=igvlog)
#END

def tepid_te_ins_igv_batch(bed, igvb, outdir):
    from collections import deque
    from subprocess import run
    indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/'
    BAMS = deque()
    M = 200
    for s in ('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15'):
        BAMS.append(indir + '{s}/{s}.bam'.format(s=s))
        # BAMS.append(indir + '{s}/{s}.discordants.sorted.bam'.format(s=s))
    GENOME = "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta"
    HEIGHT = 600
    IGV = '/srv/netscratch/dep_mercier/grp_schneeberger/software/igv/IGV_Linux_2.10.2/igv.sh'
    OUTDIR = outdir

    with open(bed, 'r') as fin, open(igvb, 'w') as fout:
        fout.write("new"+"\n")
        fout.write("genome {}\n".format(GENOME))
        for bam in BAMS:
            fout.write("load {}\n".format(bam))
        fout.write("snapshotDirectory {}\n".format(OUTDIR))
        fout.write("maxPanelHeight {}\n".format(HEIGHT))
        fout.write("squish\n")
        for line in fin:
            line = line.strip().split()
            c = line[0]
            s = min(int(line[1]), int(line[9])) - M
            e = max(int(line[2]), int(line[10])) + M
            fout.write("goto {}:{}-{}\n".format(c, s, e))
            fout.write("snapshot {}:{}-{}.png\n".format(c, s, e))
        fout.write("exit\n")
    with open(indir+"igv.log", 'a') as igvlog:
        run("{} -b {}".format(IGV, igvb).split(), stdout=igvlog, stderr=igvlog)
#END

def tepid_te_ins_bs_igv_batch(bed, igvb, outdir):
    " Plot branch specific TE ins"
    from collections import deque
    from subprocess import run
    indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/'
    BAMS = deque()
    M = 200
    for s in ('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15'):
        # BAMS.append(indir + '{s}/{s}.bam'.format(s=s))
        BAMS.append(indir + '{s}/{s}.discordants.sorted.bam'.format(s=s))
    GENOME = "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta"
    HEIGHT = 600
    IGV = '/srv/netscratch/dep_mercier/grp_schneeberger/software/igv/IGV_Linux_2.10.2/igv.sh'
    OUTDIR = outdir

    with open(bed, 'r') as fin, open(igvb, 'w') as fout:
        fout.write("new"+"\n")
        fout.write("genome {}\n".format(GENOME))
        for bam in BAMS:
            fout.write("load {}\n".format(bam))
        fout.write("snapshotDirectory {}\n".format(OUTDIR))
        fout.write("maxPanelHeight {}\n".format(HEIGHT))
        fout.write("ColorBy READ_STRAND\n")
        fout.write("squish\n")
        for line in fin:
            line = line.strip().split()
            c = line[0]
            s = int(line[1])
            e = int(line[2])
            fout.write("goto {}:{}-{}\n".format(c, s, e))
            fout.write("snapshot {}:{}-{}.png\n".format(c, s, e))
        fout.write("exit\n")
    with open(indir+"igv.log", 'a') as igvlog:
        run("{} -b {}".format(IGV, igvb).split(), stdout=igvlog, stderr=igvlog)
#END

def tefinder_igv_batch(pos, igvb, outdir):
    from collections import deque
    from subprocess import run
    indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/pheno_mut/TE_indel/tefinder/'
    BAMS = deque()
    for s in ('WT_1', 'WT_19', 'MUT_11_1', 'MUT_15'):
        BAMS.append(indir + '{s}/all_discordants.bam'.format(s=s))
    GENOME = "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta"
    HEIGHT = 600
    IGV = '/srv/netscratch/dep_mercier/grp_schneeberger/software/igv/IGV_Linux_2.10.2/igv.sh'
    OUTDIR = outdir

    with open(igvb, 'w') as fout:
        fout.write("new"+"\n")
        fout.write("genome {}\n".format(GENOME))
        for bam in BAMS:
            fout.write("load {}\n".format(bam))
        fout.write("snapshotDirectory {}\n".format(OUTDIR))
        fout.write("ColorBy READ_STRAND\n")
        fout.write("maxPanelHeight {}\n".format(HEIGHT))
        fout.write("squish\n")
        for p in pos:
            fout.write("goto {}\n".format(p))
            fout.write("snapshot {}.png\n".format(p))
        fout.write("exit\n")
    with open(indir+"igv.log", 'a') as igvlog:
        run("{} -b {}".format(IGV, igvb).split(), stdout=igvlog, stderr=igvlog)
#END

