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
    for k,v in out.items():
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

