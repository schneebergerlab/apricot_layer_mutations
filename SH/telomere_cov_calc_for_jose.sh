## SET PATHS
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/telomer_length/
telfa=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/telomeric_sequence.fa
refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta

## Align leaf re-sequencing reads to
INDIR=/biodata/dep_mercier/grp_schneeberger/reads/Apricot/Reseq-RP.P-B.P/data/seqs/
for s in A B C D; do
    cd $CWD
    mkdir leaf; cd leaf
    bsub -q ioheavy -n 10 -R "span[hosts=1] rusage[mem=25000]" -M 30000 -oo leaf_${s}.log -eo leaf_${s}.err "
    bowtie2 --end-to-end \
        --very-sensitive \
        --threads 10 \
        -x $telfa \
        -1 ${INDIR}3649_${s}_run531_*_L008_R1_001.fastq.gz \
        -2 ${INDIR}3649_${s}_run531_*_L008_R2_001.fastq.gz \
        --rg-id ${s}\"\\t\"PL:ILLUMINA\"\\t\"SM:1 \
        | samtools view -h -@20 -F4 - \
        | samtools sort -@10 -O BAM - \
        > leaf_${s}_vs_telo.bam
    samtools index -@10 leaf_${s}_vs_telo.bam
    "
    
    ### IF YOU HAVE THESE ALIGNMENTS ALREADY THEN YOU DO NOT NEED TO RE-ALIGN. OR YOU CAN REDO THE ALIGNMENT TO ENSURE CONSISTENCY
    bsub -q ioheavy -n 10 -R "span[hosts=1] rusage[mem=25000]" -M 30000 -oo leaf_${s}.log -eo leaf_${s}.err "
    bowtie2 --end-to-end \
        --very-sensitive \
        --threads 10 \
        -x $refcur \
        -1 ${INDIR}3649_${s}_run531_*_L008_R1_001.fastq.gz \
        -2 ${INDIR}3649_${s}_run531_*_L008_R2_001.fastq.gz \
        --rg-id ${s}\"\\t\"PL:ILLUMINA\"\\t\"SM:1 \
        | samtools view -h -@20 -F4 - \
        | samtools sort -@10 -O BAM - \
        > leaf_${s}_vs_cur.bam
    samtools index -@10 leaf_${s}_vs_cur.bam
    "
done

### Leaf cov count
for s in A B C D; do
    cd $CWD/leaf
    bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 10000 -oo leaf_${s}_cov.log -eo leaf_${s}_cov.err "
    syri3.8
    /srv/biodata/dep_mercier/grp_schneeberger/software/hometools/myUsefulFunctions.py bamcov leaf_${s}_vs_telo.bam leaf_${s}_vs_telo.cov &
    /srv/biodata/dep_mercier/grp_schneeberger/software/hometools/myUsefulFunctions.py bamcov leaf_${s}_vs_cur.bam leaf_${s}_vs_cur.cov
    "
done

