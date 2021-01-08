import argparse


def get_stats(D, B, P):
    import os
    from collections import deque
    from pandas import DataFrame
    import numpy as np
    with open(B, 'r') as cf:
        with open(P, 'w') as fout:
            for line in cf:
                c, cs = line.strip().split()
                print(c)
                # c = 'utg000109l'
                # cs = 1797241
                with open(D + "/" + c + ".oblast", 'r') as f:
                    s = deque()
                    e = deque()
                    gc = deque()
                    for l in f:
                        if l[0] == "#": continue
                        l = l.strip().split()
                        if int(l[3]) < 5000: continue
                        if float(l[10]) > 0.0001: continue
                        if float(l[2]) < 90: continue
                        gc.append(l[1])
                        s.append(int(l[6]))
                        e.append(int(l[7]))
                    df = DataFrame({'gc': gc,
                                    's': s,
                                    'e': e})

                    chr_ranges = {}
                    for gcid, grp in df.groupby('gc'):
                        group = grp.sort_values(['e'], ascending=False).copy()
                        group.sort_values(['s'], inplace=True)
                        ranges = deque()
                        cur_s = 0
                        cur_e = 0
                        for row in group.itertuples(index=False):
                            if row.s > cur_e:
                                ranges.append([cur_s, cur_e])
                                cur_s = row.s
                                cur_e = row.e
                            elif row.e < cur_e: continue
                            else: cur_e = row.e
                        ranges.append([cur_s, cur_e])
                        # print(ranges)
                        chr_ranges[gcid] = list(ranges)[1:]

                    keys = list(chr_ranges.keys())
                    fout.write(c + "\t" + \
                               str(cs) + "\t" + \
                               ";".join([k+":"+str(sum([v[1] - v[0] + 1 for v in chr_ranges[k]])) for k in keys]) + "\t" + \
                               ";".join([k+":"+','.join([str(v[0]) + '-' + str(v[1]) for v in chr_ranges[k]]) for k in keys]) +"\n")


def get_stats_plot(f, fout, fout2):
    from matplotlib import pyplot as plt
    from matplotlib.patches import Rectangle
    import os
    from collections import deque
    from pandas import DataFrame
    import numpy as np
    c = deque()
    cs = deque()
    mc = deque()
    mcv = deque()
    tcv = deque()
    c_reg = {}
    with open(f, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            c.append(line[0])
            cs.append(int(line[1]))
            if len(line) > 2:
                cv = {i.split(':')[0]: int(i.split(':')[1]) for i in line[2].split(';')}
                mc.append(max(cv.keys(), key=lambda x: cv[x]))
                mcv.append(int(cv[max(cv.keys(), key=lambda x: cv[x])]))
                ranges = np.sort([(int(j.split('-')[0]), int(j.split('-')[1])) for i in line[3].split(';') for j in i.split(':')[1].split(',')], 0)
                pos = np.zeros(int(line[1])+1)
                for r in ranges:
                    pos[r[0]:(r[1]+1)] = 1
                tcv.append(int(sum(pos)))
                c_reg[line[0]] = {}
                for r in line[3].split(';'):
                    ref_c = r.split(':')[0]
                    reg = [(int(i.split('-')[0]), int(i.split('-')[1])) for i in r.split(':')[1].split(',')]
                    c_reg[line[0]][ref_c] = reg
            else:
                mc.append("")
                mcv.append(0)
                tcv.append(0)

    df = DataFrame({'c'  : c,
                    'cs' : cs,
                    'mc' : mc,
                    'mcv': mcv,
                    'tcv': tcv
                    })

    df['is_bad'] = (np.array(mcv)/np.array(cs) > 0.95) | (np.array(cs) - np.array(mcv) < 10000)
    df.sort_values(['is_bad'], inplace=True, ascending=False)
    df['clrs'] = ['lightgrey' if i else 'black' for i in df.is_bad]
    df['percent_aligned'] = df.tcv*100/df.cs
    pltdf = df.loc[(df.tcv/df.cs < 0.95) & (df.cs - df.tcv > 5000) & (df.percent_aligned < 99) & (df.is_bad == False)].copy()
    pltdf.sort_values(['cs', 'percent_aligned'], inplace=True)
    fig = plt.figure(figsize=[20, 6])
    ax = fig.add_subplot(1, 1, 1)
    ax.set_axisbelow(True)
    ax.set_yticks(np.arange(0, max(df.cs), 100000))
    ax.grid(b=True, which='major', axis='y', linestyle='--')
    ax.bar(pltdf.c, pltdf.cs, label='contig_length')
    ax.bar(pltdf.c, pltdf.tcv, label='total_blast_length')
    ax.bar(pltdf.c, pltdf.mcv, label='max_blast_length')
    ax.legend()
    ax.set_xticklabels(pltdf.c, rotation=90)
    plt.subplots_adjust(left=0.02, bottom=0.18, right=0.98, top=0.98, wspace=0.0, hspace=0.25)
    fig.savefig(fout)

    from matplotlib.backends.backend_pdf import PdfPages
    pltdf = df.loc[(df.is_bad == False)].copy()
    pltdf.sort_values(['c'], inplace=True)
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
    # with PdfPages('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/TMP_blast_save_contigs.pdf') as pdf:
    with PdfPages(fout2) as pdf:
        fig = plt.figure(figsize=[10,10])
        i = 0
        for c in pltdf.c:
            # print(i)
            i += 1
            ax = fig.add_subplot(5, 1, i)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.tick_params(left=False)
            ax.set_yticklabels([])
            ax.set_xlim([0, list(pltdf.loc[pltdf.c == c].cs)[0]])
            ax.set_xlabel(c + " Non-aligning base count: {}".format(list(pltdf.loc[pltdf.c == c].cs)[0] - list(pltdf.loc[pltdf.c == c].tcv)[0]))
            ax.add_patch(Rectangle((0,0), list(pltdf.loc[pltdf.c == c].cs)[0], 1,
                                   edgecolor = 'black',
                                   fill=False,
                                   lw=3))

            color_index = 0
            try:
                for k, v in c_reg[c].items():
                    j = 0
                    for values in v:
                        if j == 0:
                            ax.add_patch(Rectangle((values[0],0 + 0.1*color_index), values[1]-values[0], 0.1,
                                                   facecolor=colors[color_index],
                                                   fill=True,
                                                   label=k,
                                                   # alpha=0.5,
                                                   lw=2))
                        else:
                            ax.add_patch(Rectangle((values[0],0 + 0.1*color_index), values[1]-values[0], 0.1,
                                                   facecolor=colors[color_index],
                                                   fill=True,
                                                   # alpha=0.5,
                                                   lw=2))
                        j+=1
                    color_index += 1
                ax.legend(loc=[1.01, 0], ncol=2)
            except KeyError as e:
                pass
            if i==5:
                plt.subplots_adjust(left=0.01, bottom=0.05, right=0.7, top=0.99, wspace=0.1, hspace=0.8)
                pdf.savefig()
                plt.close()
                fig=plt.figure(figsize=[10,10])
                i=0
                # break
        plt.subplots_adjust(left=0.01, bottom=0.05, right=0.7, top=0.99, wspace=0.1, hspace=0.8)
        pdf.savefig()
        plt.close()

def minimap2_filt(indir, fout, fout2):
    import os
    import sys
    import numpy as np
    from collections import deque
    sys.path.insert(0, '/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/')
    from myUsefulFunctions import mergeRanges

    INDIR   = indir
    FOUT    = fout
    FOUT2   = fout2
    #
    # INDIR = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/minimap2_check/'
    # FOUT  = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/minimap2_check/contigs.overlap.pdf'
    # FOUT2  = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/minimap2_check/contigs.low_overlap.txt'

    pafs = [i for i in os.listdir(INDIR) if "fasta.paf" in i]

    outdict = {}
    for paf in pafs:
        id = paf.split('.')[0]
        with open(INDIR + paf, 'r') as f:
            r = deque()
            for line in f:
                line = line.strip().split()
                if int(line[9])/int(line[10])>=0.9 and int(line[9])>1000:
                    r.append([int(line[2]), int(line[3])])
            l = int(line[1])
            r = mergeRanges(np.array(r))
            try:
                m = sum(r[:,1] - r[:, 0] + 1)
            except IndexError as e:
                m = 0
            outdict[id] = [l, m/l, r]

    from matplotlib import pyplot as plt
    ks = [k for k in sorted(outdict.keys(), key= lambda x: outdict[x][0]) if outdict[k][1]<0.9]
    fig = plt.figure(figsize = [16,8])
    plt.bar(ks, [outdict[k][0] for k in ks])
    plt.bar(ks, [outdict[k][1] * outdict[k][0]  for k in ks])
    plt.xticks(rotation=90)
    plt.xlabel("contig")
    plt.ylabel("length (total and overlapping)")
    plt.subplots_adjust(left=0.02, bottom=0.15, right=0.98, top=0.98, wspace=0.0, hspace=0.25)
    plt.savefig(FOUT)
    plt.close()

    with open(FOUT2, 'w') as f:
        for k, v in outdict.items():
            if v[1] < 0.9:
                try:
                    f.write(k + "\t" + str(v[1]) + "\t" + ";".join(["-".join(list(map(str,i))) for i in v[2]]) + "\n")
                except TypeError as e:
                    print(k, v)
                    # return


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Get stats for bad contigs blast result to good contigs")
    parser.add_argument('d', help='directory containing blast results. One blast output file per contig.', type=str)
    parser.add_argument('b', help='path to bad contigs sizes', type=argparse.FileType('r'))
    # parser.add_argument('g', help='path to good contigs sizes', type=argparse.FileType('r'))
    # parser.add_argument('-m', help='minimum read coverage cutoff', default=0, type=int)
    # parser.add_argument('-M', help='maximum read coverage cutoff', default=1000, type=int)
    # parser.add_argument('-s', help='sample IDs separated by \':\'', type=str)
    parser.add_argument('-p', dest='p', help='prefix for output file', type=str, default='')

    args = parser.parse_args()

    f = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/run3_using_p_utg/mm2_good_vs_bad_contigs.paf'

    outdict = {}
    with open(f, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            if int(line[9]) > 10000:
                if int(line[9])/int(line[10]) > 0.95:
                    try:
                        outdict[line[0]][1].append(line)

            break

