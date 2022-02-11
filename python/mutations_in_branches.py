from collections import deque, Counter
from tqdm import tqdm

## Read snps/indels between assemblies identified by syri
syriout = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.out'
syri_snp_list = {}
syri_indel_list = {}
with open(syriout, 'r') as fin:
    types = {'SNP', 'INS', 'DEL'}
    for line in fin:
        line = line.strip().split()
        if line[10] == 'SNP':
            if line[3] == 'N' or line[4] == 'N': continue
            syri_snp_list[(line[0], int(line[1]))] = (line[3], line[4])
        elif line[10] == 'INS':
            syri_indel_list[(line[0], int(line[1]))] = ('+', line[4][1:])
        elif line[10] == 'DEL':
            syri_indel_list[(line[0], int(line[1])+1)] = ('-', line[3][1:])

## Conditions for selecting variants:
## Only alt alleles with more than 5 reads are considered
## If a position has multiple alt_alleles => noise_pos_list AND not a somatic mutation
## For syri SNP/indel positions: record alt_allele_frequency AND not a somatic mutation
## If not above => save candidate SNP position
indir = '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling'
SAMPLES = ('WT_1', 'wt7', 'wt18', 'WT_19', 'mut4', 'MUT_11_1', 'mut11_2', 'MUT_15')
BASE_DICT = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
BASE_DICT2 = {4: 'A', 5: 'C', 6: 'G', 7: 'T'}
S_COV = {'WT_1': (50, 230),
         'wt7': (40, 150),
         'wt18': (40, 150),
         'WT_19': (50, 240),
         'mut4': (40, 140),
         'MUT_11_1': (50, 250),
         'mut11_2': (40, 140),
         'MUT_15': (50, 230)}
noise_pos = set()
snp_pos = set(syri_snp_list.keys())
indel_pos = set(syri_indel_list.keys())
samples_sm = {}
sample_syri_snp = {}
sample_syri_indel = {}
for sample in SAMPLES:
    sample_noise = deque()
    snp_alfrq = {}
    indel_alfrq = {}
    sm_pos = {}
    with open(f'{indir}/{sample}/filtered_low_ref_al_bam_read_counts_b30_q10.bt2.txt') as fin:
        # count = 1
        for line in tqdm(fin):
            # count += 1
            # if count == 10000:
            #     break
            line = line.strip().split()
            # Get allele frequency if position is in syri_snp
            if (line[0], int(line[1])) in snp_pos:
                p = (line[0], int(line[1]))
                snp_alfrq[p] = round(int(line[BASE_DICT[syri_snp_list[p][1]]])/int(line[3]), 2)
                continue
            # Get allele frequency if position is in syri_indel
            if (line[0], int(line[1])) in indel_pos:
                p = (line[0], int(line[1]))
                v = syri_indel_list[p][0] + syri_indel_list[p][1]
                for j in range(9, len(line), 2):
                    if line[j] == v:
                        indel_alfrq[p] = round(int(line[j+1])/int(line[3]), 2)
                        break
                continue
            # Check if the position is noisy, if not then select somatic mutation
            if (line[0], int(line[1])) not in noise_pos:
                # Check read-depth at the position
                if not S_COV[sample][0] < int(line[3]) < S_COV[sample][1]:
                    sample_noise.append((line[0], int(line[1])))
                ind = [4, 5, 6, 7] + list(range(10, len(line), 2))
                hf = 0
                for i in ind:
                    try:
                        if int(line[i]) >= 5:
                            hf += 1
                    except IndexError as e:
                        print(f'ERROR in reading line: {line}')
                if int(line[BASE_DICT[line[2]]]) >= 5:
                    hf -= 1
                if hf > 1:
                    sample_noise.append((line[0], int(line[1])))
                elif hf == 1:
                    ## Get alternate allele readcount and alternate allele ratio
                    rc = 0
                    ar = 0
                    alt = ''
                    for i in ind:
                        if int(line[i]) >= 5:
                            if BASE_DICT[line[2]] != i:
                                rc = int(line[i])
                                ar = round(int(line[i])/int(line[3]), 2)
                                try:
                                    alt = BASE_DICT2[i]
                                except KeyError:
                                    alt = line[i-1]
                    sm_pos[(line[0], int(line[1]))] = (line[2], alt, rc, ar)
    noise_pos = noise_pos.union(set(sample_noise))
    samples_sm[sample] = sm_pos
    sample_syri_snp[sample] = snp_alfrq
    sample_syri_indel[sample] = indel_alfrq

# Remove noise positions from candidate lists
for sample in SAMPLES:
    print(sample, len(samples_sm[sample]))
    for p in noise_pos:
        try:
            samples_sm[sample].pop(p)
        except KeyError:
            pass
    print(sample, len(samples_sm[sample]))

# Remove positions based on candidates
BG_SAMPLE = {'WT_1': ('WT_19', 'mut4', 'MUT_11_1', 'mut11_2', 'MUT_15'),
             'wt7': ('WT_19', 'mut4', 'MUT_11_1', 'mut11_2', 'MUT_15'),
             'wt18': ('WT_19', 'mut4', 'MUT_11_1', 'mut11_2', 'MUT_15'),
             'WT_19': ('WT_1', 'wt7', 'wt18', 'mut4', 'MUT_11_1', 'mut11_2', 'MUT_15'),
             'mut4': ('WT_1', 'wt7', 'wt18', 'WT_19'),
             'MUT_11_1': ('WT_1', 'wt7', 'wt18', 'WT_19'),
             'mut11_2': ('WT_1', 'wt7', 'wt18', 'WT_19'),
             'MUT_15': ('WT_1', 'wt7', 'wt18', 'WT_19')}

to_pop = {}
for sample in SAMPLES:
    bgpos = set()
    for bgs in BG_SAMPLE[sample]:
        bgpos = bgpos.union(set(samples_sm[bgs].keys()))
    sample_pop = deque()
    for k, v in tqdm(samples_sm[sample].items()):
        for bgs in BG_SAMPLE[sample]:
            try:
                bgv = samples_sm[bgs][k]
            except KeyError:
                continue
            if bgv[3] > 0.2*v[3]:
                sample_pop.append(k)
                break
    to_pop[sample] = sample_pop

for sample in SAMPLES:
    for pos in to_pop[sample]:
        samples_sm[sample].pop(pos)

all_pos = deque()
for sample in SAMPLES:
    all_pos.extend(list(samples_sm[sample].keys()))
all_pos = Counter(all_pos)

to_pop2 = {}
for sample in SAMPLES:



    print(sample, len(samples_sm[sample]), len(to_pop[sample]))






