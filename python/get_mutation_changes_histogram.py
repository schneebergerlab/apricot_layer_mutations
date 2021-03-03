import argparse

def count_mut_type(fin, rc, qc):
    from collections import deque
    import sys

    nuc = ('A', 'C', 'G', 'T')
    mut = {'AC': 0, 'AG': 0, 'AT': 0,
           'CA': 0, 'CG': 0, 'CT': 0,
           'GA': 0, 'GC': 0, 'GT': 0,
           'TA': 0, 'TC': 0, 'TG': 0}
    for line in fin:
        line = line.strip().split()
        try:
            if line[rc-1].upper() not in nuc or line[qc-1].upper() not in nuc:
                raise ValueError('Row {} have incorrect base. Only A, C, G, T are accepted'.format(' '.join(line)))
                sys.exit()
            mut[line[rc-1].upper() + line[qc-1].upper()] += 1
        except IndexError:
            raise IndexError('Row {} do not have columns corresponding to reference/query allele'.format(' '.join(line)))
            sys.exit()
    return mut


def get_mut_plt(args):
    if args.samples is not None:
        if len(args.samples) != len(args.f): sys.exit('Number of sample names should equal number of input files')
        SAMPLES = [s for s in args.samples]
    else:
        SAMPLES = ['Sample_' + str(i) for i in range(len(args.f))]

    if args.c is not None:
        if len(args.c) != len(args.f): sys.exit('Number of colours should equal number of input files')
        COLS = [c for c in args.c]
    else:
        from matplotlib import cm
        COLS = cm.get_cmap('Dark2').colors


    mutdata = {}
    for i in range(len(args.f)):
        fin = open(args.f[i].name, 'r')
        mutdata[SAMPLES[i]] = count_mut_type(fin, args.rc, args.qc)
        fin.close()

    from matplotlib import pyplot as plt
    import numpy as np

    fig = plt.figure(figsize=[args.W, args.H])
    ax = fig.add_subplot()
    ax.set_ylim([0, args.ymax])
    TW = 0.7        # Total width of the bars
    width = TW/len(SAMPLES)
    x = np.arange(6)
    x_off = -((len(SAMPLES)-1)/2)*width
    labels = ['A:T-->G:C',
              'G:C-->A:T',
              'A:T-->T:A',
              'G:C-->T:A',
              'A:T-->C:G',
              'G:C-->C:G']
    for i in range(len(SAMPLES)):
        y = np.array([mutdata[SAMPLES[i]]['AG'] + mutdata[SAMPLES[i]]['TC'],
                      mutdata[SAMPLES[i]]['GA'] + mutdata[SAMPLES[i]]['CT'],
                      mutdata[SAMPLES[i]]['AT'] + mutdata[SAMPLES[i]]['TA'],
                      mutdata[SAMPLES[i]]['GT'] + mutdata[SAMPLES[i]]['CA'],
                      mutdata[SAMPLES[i]]['AC'] + mutdata[SAMPLES[i]]['TG'],
                      mutdata[SAMPLES[i]]['GC'] + mutdata[SAMPLES[i]]['CG']])
        ylab = '#mutations'
        if not args.n:
            y = y*100/np.sum(y)
            ylab = 'Percentage'
        ax.bar(x + x_off, y, width, label=SAMPLES[i], color=COLS[i])
        x_off += width

    ax.legend()
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel(ylab)
    ax.set_title('Mutation changes')
    plt.tight_layout()
    plt.savefig(args.o.name)
    plt.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Generate histograms for different mutation changes using table input", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", help="Input file containing the frequency data", type=argparse.FileType('r'), nargs='+')
    parser.add_argument("-rc", help='Column index containing reference allele', type=int, default=1)
    parser.add_argument("-qc", help='Column index containing query allele', type=int, default=2)
    parser.add_argument("-o", help="Output file name", type=argparse.FileType('w'), default='mutant_change_hist.pdf')
    parser.add_argument("-samples", help='Sample names to use in plot', type=str, nargs='*')
    parser.add_argument("-n", help='Y-axis corresponds to mutations counts (and not percentages)', action='store_true', default=False)
    parser.add_argument("-W", help="width of the plot (in inches)", type=float, default=6)
    parser.add_argument("-H", help="height of the plot (in inches)", type=float, default=4)
    parser.add_argument("-c", help='color associated with different sample', type=str, nargs='*')
    parser.add_argument("-t", help="title of the plot", type=str, default='Mutation changes', nargs='+')
    parser.add_argument("-x", help="X-axis label", type=str, default='value', nargs='+')
    parser.add_argument("-y", help="Y-axis label", type=str, default='counts', nargs='+')
    parser.add_argument("-ymax", help="Upper limit for y-axis", type=float, default=None)
    # parser.add_argument("-xlog", help="X-axis on log scale", default=False, action='store_true')
    # parser.add_argument("-ylog", help="Y-axis on log scale", default=False, action='store_true')


    args = parser.parse_args()
    # args = parser.parse_args('-f /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/MUT_11_1/TMP1.txt /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/MUT_11_1/TMP2.txt /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/MUT_11_1/TMP2.txt -rc 6 -qc 7 -samples TMP1 TMP2 TMP2 -ymax 40'.split())

    get_mut_plt(args)