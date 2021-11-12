#!/usr/bin/env python3
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser("Creates a batch file that can be used with IGV to get snapshot for specific regions", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("pos", help='Position', type=int)
    parser.add_argument("chr", help='Chromosome', type=str)
    parser.add_argument("-t", help="Plot DNA or RNA bam", type=str, choices=['dna', 'rna'], default='dna')
    parser.add_argument("-m", help='Margin around position', type=int, default=1000)
    parser.add_argument("-p", help='Number of vertical pixels', type=int, default=300)
    parser.add_argument("-o", help="Output BATCH file name", type=argparse.FileType('w'), default='igv.bat')
    parser.add_argument("-d", help="Output directory for snapshots", type=str, default='.')
    parser.add_argument("--run", help="Run IGV batch command", default=False, action="store_true")
    parser.add_argument("--igv", help="Path to igv.sh", type=argparse.FileType('r'))


    args = parser.parse_args()

    import os
    if not os.path.isdir(args.d):
        try:
            os.mkdir(args.d)
        except:
            raise argparse.ArgumentTypeError("{} is not a valid path".format(args.d))

    SAMPLES = ['WT_1', 'WT_19', 'MUT_11_1', 'MUT_15']
    if args.t == 'dna':
        BAMS = ['/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/{s}/{s}.sorted.bt2.bam'.format(s=s) for s in SAMPLES]
    elif args.t == 'rna':
        BAMS = ['/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/get_cells/{s}/{s}/outs/possorted_genome_bam.bam'.format(s=s) for s in SAMPLES]
    else:
        raise argparse.ArgumentError("Incorrect bam type")
    POS = args.pos
    CHR = args.chr
    WIDTH = args.m
    HEIGHT = args.p
    OUTFILE = args.o.name
    OUTDIR = args.d
    GENOME = "/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta"
    RUN = args.run
    IGV = args.igv.name if args.igv is not None else None

    with open(OUTFILE, 'w') as fout:
        fout.write("new"+"\n")
        fout.write("genome {}\n".format(GENOME))
        for bam in BAMS:
            fout.write("load {}\n".format(bam))
        fout.write("snapshotDirectory {}\n".format(OUTDIR))
        fout.write("maxPanelHeight {}\n".format(HEIGHT))
        fout.write("squish\n")
        fout.write("goto {}:{}-{}\n".format(CHR, POS-WIDTH, POS+WIDTH))
        fout.write("group BASE_AT_POS {}:{}\n".format(CHR, POS))
        fout.write("snapshot {}:{}-{}.png\n".format(CHR, POS-WIDTH, POS+WIDTH))
        fout.write("exit\n")
    if RUN:
        if IGV is not None:
            igv = IGV
        else:
            igv = '/srv/netscratch/dep_mercier/grp_schneeberger/software/igv/IGV_Linux_2.10.2/igv.sh'
        from subprocess import run
        with open("igv.log", 'a') as igvlog:
            run("{} -b {}".format(igv, OUTFILE).split(), stdout=igvlog, stderr=igvlog)
