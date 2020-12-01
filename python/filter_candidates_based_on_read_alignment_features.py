import argparse
from collections import deque
QS= "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~"


if __name__ == '__main__':
    """
    For each candidate sSNV position, check the mpileup data to see whether is there:
        1) Read strand bias, i.e. reads on the position are biased towards one strand (Using binomial test)
        2) Allele strand bias, i.e. reference and alternate alleles have different distribution on the strands (Using binomial test)
        3) Base quality of the alternate alleles different from base quality of the reference alleles (Using Mann–Whitney U test)
        4) Mapping quality of the alternate alleles different from mapping quality of the reference alleles (Using Mann–Whitney U test)
        5) Alignment Score of the reads with alternate alleles too different from the reads with reference alleles (Using Mann–Whitney U test)
        6) Read position bias, i.e. the variant bases are on one position on a read
    
    Check the SAM file to see:
        1) If the positions is represented by both Read1 and Read2
        
    
    """
    parser = argparse.ArgumentParser("Select positions which are supported by read-level features",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('reg', help='List of regions to check', type=argparse.FileType('r'))
    parser.add_argument('bam', help='bam file', type=argparse.FileType('r'))
    parser.add_argument('ref', help='path to reference file', type=argparse.FileType('r'))

    parser.add_argument('-n', dest='n', help='minimum alternate read count to select position', type=int, default=1)
    parser.add_argument('-o', dest='o', help='prefix for the output file', default='cell_var', type=str)

    args = parser.parse_args()

    REG = args.reg.name
    BAM = args.bam.name
    OUT = args.
    q

    REG = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_1/multi_cell_WT_1_only_SNPs_candidate.sorted.common.unique_genome.non_repeat.non_indel.regions'
    BAM = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_1/WT_1.sorted.RG.bam'
    REF = '/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
    OUT = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_1/test.pileup'

    from subprocess import Popen, PIPE

    # Get MPILEUP DATA
    cmd = '/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_manish/samtools mpileup -A -q 30 -Q 10 -E --reference {PATH_TO_REF} -r {REGION} --ff UNMAP,SECONDARY,QCFAIL,DUP -s --output-QNAME --output-extra AS {PATH_TO_BAM}'
    print(cmd)
    count = 0
    variants = deque()
    with open(REG, 'r') as fin:
        # with open(OUT, 'w') as fout:
        for line in fin:
            # count += 1
            # if count == 10 : break
            line = line.strip().split()
            variants.append(line)
                region = '{}:{}-{}'.format(line[0], line[1], line[2])
                run = cmd.format(PATH_TO_REF=REF, REGION=region, PATH_TO_BAM=BAM)
                p = Popen( run.split(), stdout=PIPE, stderr=PIPE)
                fout.write(p.communicate()[0].decode())

    variants = tuple(variants)

    class snvdata:
        def __init__(self, ls):
            self.chr = ls[0]
            self.pos = ls[1]
            self.ref = ls[2]
            self.rc  = int(ls[3])
            self.indelcount, self.bases = self._getbases(ls[4])
            if len(self.bases) != self.rc:
                raise Exception('Number of identified bases if not equals to read count for {}:{}'.format(self.chr, self.pos))
            self.BAQ = [ord(c) - 33 for c in ls[5]]
            self.MQ = [ord(c) - 33 for c in ls[6]]
            self.AS = list(map(int, ls[8].split(',')))


        def _getbases(self, l):
            bases = deque()
            skip = 0
            indel = False
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

        def getBAQ(self, base):
            return [self.BAQ[i] for i in range(len(self.bases)) if self.bases[i] == base]

        def getMQ(self, base):
            return [self.MQ[i] for i in range(len(self.bases)) if self.bases[i] == base]

        def getAS(self, base):
            return [self.AS[i] for i in range(len(self.bases)) if self.bases[i] == base]



    snvs = deque()
    ## Read data from the mpileup file
    count = 0
    with open(OUT, 'r') as fin:
        for line in fin:
            line = line.strip().split()
            snvs.append(snvdata(line))


    from scipy.stats import binom_test, fisher_exact, mannwhitneyu

    stats = deque()
    count = 0
    for i in range(len(variants)):
        # Count alternate allele count on each strand
        upcnt = snvs[i].basecnt(variants[i][6].upper())
        locnt = snvs[i].basecnt(variants[i][6].lower())

        # # Do binomial test on strand
        pbin = binom_test(upcnt, upcnt + locnt, 0.5, 'two-sided')

        # Perform fisher's exact test to test for strand bias
        rf = snvs[i].basecnt('.')
        rr = snvs[i].basecnt(',')
        af = upcnt
        ar = locnt
        pfis = fisher_exact([[rf, af], [rr, ar]])[1]

        # Do Mann-Whitney U test for BAQ
        a = snvs[i].getBAQ(variants[i][6].upper()) + snvs[i].getBAQ(variants[i][6].lower())
        b = snvs[i].getBAQ('.') + snvs[i].getBAQ(',')
        pbak = mannwhitneyu(a, b, alternative='two-sided').pvalue

        # Do Mann-Whitney U test for MQ
        a = snvs[i].getMQ(variants[i][6].upper()) + snvs[i].getMQ(variants[i][6].lower())
        b = snvs[i].getMQ('.') + snvs[i].getMQ(',')
        try:
            pmq = mannwhitneyu(a, b, alternative='two-sided').pvalue
        except ValueError as e:
            # x = e
            # print(type(e), e.message)
            if e.__str__() == 'All numbers are identical in mannwhitneyu':
                pmq = 1
            else:
                raise Exception(e)

        # Do Mann-Whitney U test for AS
        a = snvs[i].getAS(variants[i][6].upper()) + snvs[i].getAS(variants[i][6].lower())
        b = snvs[i].getAS('.') + snvs[i].getAS(',')
        pas = mannwhitneyu(a, b, alternative='two-sided').pvalue


        stats.append([pbin, upcnt, locnt, pfis, pbak, pmq, pas])

    # pvbin = np.array([True if j < 0.05 else False for j in p_adjust([i[0] for i in stats], 'BH')])
    pvbin = np.array([True if j < 0.05 else False for j in [i[0] for i in stats]])
    notup = np.array([True if j == 0   else False for j in [i[1] for i in stats]])
    notlo = np.array([True if j == 0   else False for j in [i[2] for i in stats]])
    # pvfis = np.array([True if j < 0.05 else False for j in p_adjust([i[3] for i in stats], 'BH')])
    pvfis = np.array([True if j < 0.05 else False for j in [i[3] for i in stats]])

    bad_pos = pvbin | notup | notlo | pvfis
    #
    # only_pvbin = pvbin & ~notup & ~notlo

    good_variants = [variants[i] for i in range(len(variants)) if not bad_pos[i]]

import matplotlib.pyplot  as plt
import numpy as np
import sys
sys.path.insert(1, '/srv/biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python')
from myUsefulFunctions import p_adjust

plt.yscale('log')
plt.hist([int(i[3]) for i in variants], bins=50)
plt.hist([int(i[3]) for i in good_variants], bins = 50)
plt.hist([float(i[4]) for i in good_variants])

r = (0.7341622, 0.2853933, 0.2425526, 0.1627632, 0.1764481, 0.1970707, 0.8895130, 0.1338674, 0.4563071, 0.1368662)
p_adjust(r, 'BH')
 [1] 0.8157358 0.4077048 0.4042544 0.3941414 0.3941414 0.3941414 0.8895130 0.3941414 0.5703839 0.3941414

def bh_adjust(p):
    " perform benjamini-hochberg correction on a list of p-values"

samtools mpileup -q 20 -Q 30 -O -r CUR1G:1046502-1046502 -s --output-QNAME --output-extra AS --reference /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/variantCorrected_AND_kmer_masked_d100-kmer80bp-8_4_Currot_on_Original_v1.0.fasta WT_1.sorted.RG.bt2.bam

l = '...,..,,..,.,,,.,..$.$.$.,..,,,.,..,,...,...,,...,.,...,...,..,,,.,,,,+1a...,,,....,,........,.,..,,,,,,,,....,,...,     13 ....,,,,......,,,..,...,..,.,,...,.,CCCC.^].^].'

def getbases(l):
    bases = deque()
    skip = 0
    indel = False
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