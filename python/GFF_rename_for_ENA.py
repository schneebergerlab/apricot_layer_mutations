# The GFF file used for analysis was not getting validated. Here, I filter/clean the GFF so that it can be uploaded to ENA.

# Things to fix: Overlapping (duplicated) UTRs

# NOTE: This script is setup to run on LMU servers. The final submission was done using local computer.

import os
from collections import deque
from gzip import open as gzopen
from tqdm import tqdm

# Defaults
cwd = '/home/ra98jam/projects/apricots/data'
report = 'ora.genome.v1.embl.gz.report'
embfin = 'ora.genome.v1.embl.gz'
gfffin = 'ora.pasa_out.3utr.sort.no_source.no_dup.gff3'
gffout = 'ora.pasa_out.3utr.sort.no_source.no_dup.clean.gff3'


class gffelement():
    def __init__(self, line):
        # line = line.strip().split('\t')
        self.chrom = line[0]
        self.type = line[2]
        self.start = int(line[3])
        self.end = int(line[4])
        self.strand = line[6]
        self.phase = line[7]
        self.rawtags = line[8]
        self.tags = {i.split("=")[0]: i.split("=")[1] for i in line[8].split(";")}

    def tostring(self):
        return "\t".join([self.chrom, '.', self.type, str(self.start), str(self.end), '.', self.strand, str(self.phase), self.rawtags])


class gffexon(gffelement):
    def __init__(self, line, i):
        super().__init__(line)
        self.id = i
        self.cds = None
        self.utr5 = None
        self.utr3 = None
    def addele(self, ele):
        if ele.type == 'CDS':
            self.cds = ele
        elif ele.type == 'five_prime_UTR':
            self.utr5 = ele
        elif ele.type == 'three_prime_UTR':
            self.utr3 = ele

class gffmra(gffelement):
    def __init__(self, line, i):
        super().__init__(line)
        self.id = i
        self.exons = deque()

    def addexon(self, line, i):
        self.exons.append(gffexon(line, i))

    def addele(self, line, i):
        ele = gffelement(line)
        # count number of matching exons. give error if more than one matching exons are present
        # stex = [ex for ex in self.exons if ele.start == ex.start]
        # endex = [ex for ex in self.exons if ele.end == ex.end]
        # allex = stex + endex
        allex = [ex for ex in self.exons if ex.start <= ele.start <= ex.end and ex.start <= ele.end <= ex.end]
        if len(allex) == 1:
            allex[0].addele(ele)
        else:
            raise ValueError(line)
        # else:
        #     try:
        #         assert len(set([ex.id for ex in stex] + [ex.id for ex in endex])) == 1
        #         allex[0].addele(ele)
        #     except AssertionError as e:
        #         raise AssertionError(line)
        #         # print(line, e)

class gffgene(gffelement):
    def __init__(self, line, i):
        super().__init__(line)
        self.id = i
        self.mrnas = deque()

    def addmrna(self, line, i):
        self.mrnas.append(gffmra(line, i))



# Read list of lines that are failing the validation:
problemanno = deque()
with open(f'{cwd}/{report}', 'r') as fin:
    for line in fin:
        line = line.strip().split(' ')
        problemanno.append([line[1].replace('"', ''), int(line[12].split('-')[0]), int(line[17].split('-')[0])])
problemlines = set([j for i in problemanno for j in i[1:]])

# Get GFF coordinates at the lines
lineanno = deque()
with gzopen(f'{cwd}/{embfin}', 'r') as fin:
    for i, line in enumerate(fin):
        if i+1 in problemlines:
            line = line.decode().strip().split()[1:]
            # print(line)
            if 'UTR' in line[0]:
                line[1] = line[1].replace("complement(", "").replace(")", "")
                # print(line)
                lineanno.append(tuple([line[0]] + line[1].split("..")))

lineanno = [[i[0], int(i[1]), int(i[2])] for i in set(lineanno)]
# Check if same line (UTR) is present twice
assert len(lineanno)*2 == len(set([j for i in lineanno for j in i[1:]]))

# Read all gene models
with open(f'{cwd}/{gfffin}', 'r') as fin:
    geneid = ''
    strand = ''
    genes = deque()
    curgene = ''
    curmrna = ''
    for i, line in enumerate(fin):
        if line[0] == '#': continue

        line = line.strip().split('\t')
        if line[2] == 'gene':
            if curgene == '':
                curgene = gffgene(line, i)
            else:
                genes.append(curgene)
                curgene = gffgene(line, i)
                # curmrna = ''
        elif line[2] in ['mRNA', 'RNA']:
            # if curmra == '':
            curgene.addmrna(line, i)
            curmrna = curgene.mrnas[-1]
            # else:
        elif line[2] == 'exon':
            curmrna.addexon(line, i)
        elif line[2] in ['CDS', 'five_prime_UTR', 'three_prime_UTR']:
            curmrna.addele(line, i)
    genes.append(curgene)
gcnt = deque()

for l in tqdm(lineanno):
    ann = 'utr5' if '5' in l[0] else 'utr3'
    fixgenes = deque()
    for g in genes:
        for i, m in enumerate(g.mrnas):
            for e in m.exons:
                if getattr(e, ann) is not None:
                    if getattr(e, ann).start == l[1]:
                        if getattr(e, ann).end == l[2]:
                            fixgenes.append([g, i])
    selectgene, mid = fixgenes[0]
    for g in list(fixgenes)[1:]:
        if len(selectgene.mrnas) > len(g[0].mrnas):
            selectgene, mid = g
    if len(selectgene.mrnas) == 2:
        print(mid, vars(selectgene.mrnas[0]), vars(selectgene.mrnas[1]))
    gcnt.append(len(selectgene.mrnas))
    # Manual check showed that just changing the selected mRNA would work for currot genes
    anndir = selectgene.strand

    for e in selectgene.mrnas[mid].exons:
        if getattr(e, ann) is not None:
            if ann == 'utr5' and anndir == '+':
                e.utr5.start += 1
                e.start += 1
                selectgene.mrnas[mid].start += 1
                if len(selectgene.mrnas) == 1:
                    selectgene.start += 1
                else:
                    print(f'manually check {selectgene.id}, {selectgene.tags}')
            elif ann == 'utr5' and anndir == '-':
                e.utr5.end -= 1
                e.end -= 1
                selectgene.mrnas[mid].end -= 1
                if len(selectgene.mrnas) == 1:
                    selectgene.end -= 1
                else:
                    print(f'manually check {selectgene.id}, {selectgene.tags}')
            elif ann == 'utr3' and anndir == '+':
                e.utr3.end -= 1
                e.end -= 1
                selectgene.mrnas[mid].end -= 1
                if len(selectgene.mrnas) == 1:
                    selectgene.end -= 1
                else:
                    print(f'manually check {selectgene.id}, {selectgene.tags}')
            elif ann == 'utr3' and anndir == '-':
                e.utr3.start += 1
                e.start += 1
                selectgene.mrnas[mid].start += 1
                if len(selectgene.mrnas) == 1:
                    selectgene.start += 1
                else:
                    print(f'manually check {selectgene.id}, {selectgene.tags}')


# Write the updated GFF
with open(f'{cwd}/{gffout}', 'w') as fout:
    for gene in genes:
        fout.write(gene.tostring() + "\n")
        for m in gene.mrnas:
            fout.write(m.tostring() + "\n")
            for e in m.exons:
                fout.write(e.tostring() + "\n")
                for a in ['cds', 'utr5', 'utr3']:
                    if getattr(e, a) is not None:
                        fout.write(getattr(e, a).tostring() + "\n")