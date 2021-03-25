import argparse
from collections import deque
QS= "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"


if __name__ == '__main__':
    """
    # For each candidate sSNV position, check the mpileup data to see whether is there:
    # 1) Read strand bias, i.e. reads on the position are biased towards one strand (Using binomial test)
    # 2) Allele strand bias, i.e. reference and alternate alleles have different distribution on the strands (Using fischer's exact test)
    # 3) Base quality of the alternate alleles different from base quality of the reference alleles (Using MannWhitney U test)
    # 4) Mapping quality of the alternate alleles different from mapping quality of the reference alleles (Using MannWhitney U test)
    # 5) Alignment Score of the reads with alternate alleles too different from the reads with reference alleles (Using MannWhitney U test)
    # 6) Read position bias, i.e. the variant bases are on one position on a read (using KS-test)
    # 7) If the positions is represented by both Read1 and Read2 (using binomial test)
    """


    parser = argparse.ArgumentParser("Select positions which are supported by read-level features",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('bed', help='BED file containing SNP positions', type=argparse.FileType('r'))
    parser.add_argument('bam', help='bam file', type=argparse.FileType('r'))
    parser.add_argument('ref', help='path to reference file', type=argparse.FileType('r'))
    parser.add_argument('-o', dest='o', help='prefix for the output file', type=str)

    args = parser.parse_args()
    BED = args.bed.name
    BAM = args.bam.name
    REF = args.ref.name
    OUT = '' if args.o is None else args.o + '_'


    # REG = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_1/multi_cell_WT_1_filtered_SNPs_candidate.sorted.unique_genome.non_repeat.bed'
    #
    # BAM = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_1/WT_1.sorted.RG.bam'
    # REF = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
    # OUT = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_1/test.pileup'

    from subprocess import Popen, PIPE
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    # Get MPILEUP DATA
    cmd = '/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools mpileup -A -q 10 -Q 30 -E --reference {PATH_TO_REF} -r {REGION} --ff UNMAP,SECONDARY,QCFAIL,DUP  --output-BP --output-MQ --output-extra FLAG --output-extra AS {PATH_TO_BAM}'
    # print(cmd)
    # count = 0


    # with open(OUT + 'snps.pileup', 'w') as fout:
    #     with open(BED, 'r') as fin:
    #         for line in fin:
    #             # count += 1
    #             # if count == 20 : break
    #             line = line.strip().split()
    #             # variants.append(line)
    #             region = '{}:{}-{}'.format(line[0], int(line[1])+1, line[2])
    #             run = cmd.format(PATH_TO_REF=REF, REGION=region, PATH_TO_BAM=BAM)
    #             p = Popen( run.split(), stdout=PIPE, stderr=PIPE)
    #             fout.write(p.communicate()[0].decode())
    variants = deque()
    with open(BED, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            variants.append(line)
    variants = tuple(variants)


    class snvdata:
        def __init__(self, ls):
            self.chr = ls[0]
            self.pos = ls[1]
            self.ref = ls[2]
            self.rc  = int(ls[3])
            self.indelcount, self.bases = self._getbases(ls[4])
            if len(self.bases) != self.rc:
                raise Exception('Number of identified bases if not equals to read count for {}:{}. ReadCount: {}, BaseCount: {}'.format(self.chr, self.pos, self.rc, len(self.bases)))
            self.BQ  = [ord(c) - 33 for c in ls[5]]
            self.MQ  = [ord(c) - 33 for c in ls[6]]
            self.RP  = [int(c) for c in ls[7].split(',')]
            self.F   = [int(c) for c in ls[8].split(',')]
            self.AS  = [int(c) for c in ls[9].split(',')]


        def _getbases(self, l):
            indelcount = 0
            bases = deque()
            skip = 0
            indel = False
            for c in l:
                if skip > 0 and indel == False:
                    skip -= 1
                    continue
                if indel == True:
                    if c in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']:
                        skip = (skip*10) + int(c)
                        continue
                    # skip = int(c)
                    else:
                        indel = False
                        skip -= 1
                        continue
                if c == '*':
                    # self.indelcount += 1
                    indelcount += 1
                    bases.append(c)
                    continue
                if c == '$': continue
                if c in ['<', '>']: sys.exit('spliced alignment found')
                if c == '^':
                    skip = 1
                    continue
                if c in ['+', '-']:
                    indel = True
                    continue
                bases.append(c)
            return [indelcount, list(bases)]


        def forwardcnt(self):
            try:
                return self.fcnt
            except AttributeError as e:
                self.fcnt = len([1 for i in self.bases if i in ['.', 'A', 'C', 'G', 'T']])
                self.rcnt = self.rc - self.fcnt
                return self.fcnt

        def reversecnt(self):
            try:
                return self.rcnt
            except AttributeError as e:
                self.rcnt = len([1 for i in self.bases if i in [',', 'a', 'c', 'g', 't']])
                self.fcnt = self.rc - self.rcnt
                return self.rcnt

        def basecnt(self, base):
            return len([1 for i in self.bases if i == base])

        def getBQ(self, base):
            return [self.BQ[i] for i in range(len(self.bases)) if self.bases[i] == base]

        def getMQ(self, base):
            return [self.MQ[i] for i in range(len(self.bases)) if self.bases[i] == base]

        def getRP(self, base):
            return [self.RP[i] for i in range(len(self.bases)) if self.bases[i] == base]

        def getF(self, base):
            return [self.F[i] for i in range(len(self.bases)) if self.bases[i] == base]

        def getAS(self, base):
            return [self.AS[i] for i in range(len(self.bases)) if self.bases[i] == base]


    snvs = deque()
    ## Read data from the mpileup file
    count = 0
    with open(OUT + 'snps.pileup', 'r') as fin:
        for line in fin:
            line = line.strip().split()
            snvs.append(snvdata(line))


    from scipy.stats import binom_test, fisher_exact, mannwhitneyu, kstest
    stats = deque()
    count = 0
    for i in range(len(snvs)):
        # Count alternate allele count on each strand
        upcnt = snvs[i].basecnt(variants[i][4].upper())
        locnt = snvs[i].basecnt(variants[i][4].lower())

        # # Do binomial test on strand
        pbin = binom_test(upcnt, upcnt + locnt, 0.5, 'two-sided')

        # Perform fisher's exact test to test for strand bias
        rf = snvs[i].basecnt('.')
        rr = snvs[i].basecnt(',')
        af = upcnt
        ar = locnt
        pfis = fisher_exact([[rf, af], [rr, ar]])[1]

        # Do Mann-Whitney U test for BAQ
        a = snvs[i].getBQ(variants[i][4].upper()) + snvs[i].getBQ(variants[i][4].lower())
        b = snvs[i].getBQ('.') + snvs[i].getBQ(',')
        try:
            pbak = mannwhitneyu(a, b, alternative='two-sided').pvalue
        except ValueError as e:
            if e.__str__() == 'All numbers are identical in mannwhitneyu':
                pbak = 1
            else:
                raise Exception(e)

        # Do Mann-Whitney U test for MQ
        a = snvs[i].getMQ(variants[i][4].upper()) + snvs[i].getMQ(variants[i][4].lower())
        b = snvs[i].getMQ('.') + snvs[i].getMQ(',')
        try:
            pmq = mannwhitneyu(a, b, alternative='two-sided').pvalue
        except ValueError as e:
            if e.__str__() == 'All numbers are identical in mannwhitneyu':
                pmq = 1
            else:
                raise Exception(e)

        # Do KS-test to check whether the mutant alleles have biased read position
        a = snvs[i].getRP(variants[i][4].upper()) + snvs[i].getRP(variants[i][4].lower())
        a_norm = [j/150 for j in a]
        if len(a_norm) == 0:
            pks = 0
        else:
            pks = kstest(a_norm, 'uniform').pvalue

        # Check if alternate is present in both Read1 and Read2
        flags = snvs[i].getF(variants[i][4].upper()) + snvs[i].getF(variants[i][4].lower())
        r1, r2 = 0, 0
        for flag in flags:
            if format(flag, '012b')[5] == '1': r1 += 1
            else: r2 += 1
        pf = binom_test(r1, r1 + r2, 0.5, 'two-sided')

        # # Do Mann-Whitney U test for AS
        # a = snvs[i].getAS(variants[i][4].upper()) + snvs[i].getAS(variants[i][4].lower())
        # b = snvs[i].getAS('.') + snvs[i].getAS(',')
        # pas = mannwhitneyu(a, b, alternative='two-sided').pvalue
        # stats.append([pbin, upcnt, locnt, pfis, pbak, pmq, pks, r1, r2, pf, pas])
        stats.append([pbin, upcnt, locnt, pfis, pbak, pmq, pks, r1, r2, pf])

    # pvbin = np.array([True if j < 0.05 else False for j in p_adjust([i[0] for i in stats], 'BH')])
    pvbin = np.array([True if j < 0.05 else False for j in [i[0] for i in stats]])
    notup = np.array([True if j < 1    else False for j in [i[1] for i in stats]])
    notlo = np.array([True if j < 1    else False for j in [i[2] for i in stats]])
    lessr = np.array([True if i[1] + i[2] < 5 else False for i in stats])
    # pvfis = np.array([True if j < 0.05 else False for j in p_adjust([i[3] for i in stats], 'BH')])
    pvfis = np.array([True if j < 0.05 else False for j in [i[3] for i in stats]])
    pvbak = np.array([True if j < 0.05 else False for j in [i[4] for i in stats]])
    pvmq  = np.array([True if j < 0.05 else False for j in [i[5] for i in stats]])
    pvks  = np.array([True if j < 0.05 else False for j in [i[6] for i in stats]])
    notr1 = np.array([True if j < 1    else False for j in [i[7] for i in stats]])
    notr2 = np.array([True if j < 1    else False for j in [i[8] for i in stats]])
    pvf   = np.array([True if j < 0.05 else False for j in [i[9] for i in stats]])
    # pvas  = np.array([True if j < 0.05 else False for j in [i[6] for i in stats]])

    df = pd.DataFrame({'pvbin': pvbin.astype('int'),
                       'lessr': lessr.astype('int'),
                       'pvfis': pvfis.astype('int'),
                       'pvbak': pvbak.astype('int'),
                       'pvmq': pvmq.astype('int'),
                       'pvks': pvks.astype('int'),
                       'pvf': pvf.astype('int')})
    df['fail'] = pvbin.astype('int') + lessr + pvfis + pvbak + pvmq + pvks + pvf
    df['index'] = range(df.shape[0])
    df.sort_values(['pvbin', 'lessr', 'pvfis', 'pvbak', 'pvmq', 'pvks', 'pvf'], inplace=True)
    df.index = range(df.shape[0])

    fig, ax = plt.subplots(figsize=[10, 10])
    sns.heatmap(df[['pvbin', 'lessr', 'pvfis', 'pvbak', 'pvmq', 'pvks', 'pvf', 'fail']], cmap = 'Blues', ax = ax)
    ax.set_xticklabels(['read-strand (binomial)', 'low alt count', 'read-strand (fisher)', 'Base-alignment', 'Mapping Quality', 'Position on read', 'Read1 vs Read2', 'fail'], rotation=45)
    plt.tight_layout()
    plt.savefig(OUT + 'variants_test.pdf')
    plt.close()

    df.to_csv(OUT + 'variants_test.csv', sep='\t', index=False)
    with open(OUT + 'good_variants.txt', 'w') as fout:
        for i in range(df.shape[0]):
            if df.iloc[i]['fail'] == 0:
                fout.write('\t'.join(variants[df.iloc[i]['index']]) + '\n')

