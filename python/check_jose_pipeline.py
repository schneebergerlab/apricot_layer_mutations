GENOME='/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Triobin_scratch/analysis/seqs/q_read/2018.10.19_FreqKmers/Data/5.Polishing/2018.10.24.HighSensitivity_polished/2.pilon_polished_assembly/A_Currot/pilon/A_Currot_pilon.fasta'
bamA='/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Reseq-RP.P-B.P_scratch/analysis/seqs/quality_reads-2018-09-11/RojoPasion/bowtie2_RP/indexed_A_Currot/A_Currot_pilon_bowtie2samtools_3649_4024/p4024_A/bt2_4024_A_PE.sorted.bam'
bamB='/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Reseq-RP.P-B.P_scratch/analysis/seqs/quality_reads-2018-09-11/RojoPasion/bowtie2_RP/indexed_A_Currot/A_Currot_pilon_bowtie2samtools_3649_4024/p4024_B/bt2_4024_B_PE.sorted.bam'
bamC='/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Reseq-RP.P-B.P_scratch/analysis/seqs/quality_reads-2018-09-11/RojoPasion/bowtie2_RP/indexed_A_Currot/A_Currot_pilon_bowtie2samtools_3649_4024/p4024_C/bt2_4024_C_PE.sorted.bam'
bamD='/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Reseq-RP.P-B.P_scratch/analysis/seqs/quality_reads-2018-09-11/RojoPasion/bowtie2_RP/indexed_A_Currot/A_Currot_pilon_bowtie2samtools_3649_4024/p4024_D/bt2_4024_D_PE.sorted.bam'
vars='/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Apricot/Reseq-RP.P-B.P_scratch/analysis/seqs/quality_reads-2018-09-11/RojoPasion/bowtie2_RP/indexed_A_Currot/A_Currot_pilon_bowtie2samtools_3649_4024/toReport_SHORE_4samples_qv/diMutaFilter_v2/sample_mQ20_mC10_iaf0.1/sample_mQ20_mC10_iaf0.1_grp1_specific_mutations.txt'
CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/tests/jose_out_check'
BAMS = (bamA, bamB, bamC, bamD)

import os
os.chdir(CWD)
with open(vars, 'r') as fin:
    WIDTH = 200
    HEIGHT = 300
    OUTFILE = "igv.bat"
    OUTDIR = CWD
    IGV='/srv/netscratch/dep_mercier/grp_schneeberger/software/igv/IGV_Linux_2.10.2/igv.sh'
    with open(OUTFILE, 'w') as fout:
        fout.write("new"+"\n")
        fout.write("genome {}\n".format(GENOME))
        for bam in BAMS:
            fout.write("load {}\n".format(bam))
        fout.write("snapshotDirectory {}\n".format(OUTDIR))
        fout.write("maxPanelHeight {}\n".format(HEIGHT))
        fout.write("squish\n")
        for line in fin:
            line = line.strip().split()
            CHR = line[0].split('_')[0] + "|arrow|pilon"
            POS = int(line[1])
            fout.write("goto {}:{}-{}\n".format(CHR, POS-WIDTH, POS+WIDTH))
            fout.write("group BASE_AT_POS {}:{}\n".format(CHR, POS))
            fout.write("snapshot {}:{}-{}.png\n".format(CHR, POS-WIDTH, POS+WIDTH))
        fout.write("exit\n")
with open("igv.log", 'a') as igvlog:
    run("{} -b {}".format(IGV, OUTFILE).split(), stdout=igvlog, stderr=igvlog)