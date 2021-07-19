# Link data to working folder
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_isoseq
ln -s /srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4876/4876_A_run545_CCS.bam WT_1.iso_seq.ccs.bam
ln -s /srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4876/4876_B_run545_CCS.bam WT_19.iso_seq.ccs.bam
ln -s /srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4876/4876_C_run545_CCS.bam MUT_11_1.iso_seq.ccs.bam
ln -s /srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4954/4954_A_run545_CCS.bam MUT_15.iso_seq.ccs.bam

# Select reads having good primer orientation
SAMPLES=("MUT_11_1" "MUT_15" "WT_1" "WT_19")
for sample in ${SAMPLES[@]}; do
    lima --isoseq --dump-clips --dump-removed ${sample}.iso_seq.ccs.bam 10x_primers.fasta ${sample}.iso_seq.bam
done

# Map reads
INDIR='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_isoseq/'
CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/'
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'


for sample in ${SAMPLES[@]}; do
    cd $CWD
    mkdir $sample; cd $sample
    bsub -q multicore20 -n 20 -R "span[hosts=1] rusage[mem=20000]" -M 25000 -oo align_cells.log -eo align_cells.err "
        # Get Reads
#        samtools fastq -@ 10 ${INDIR}${sample}.iso_seq.5p--3p.bam | pigz -p 10 > ${sample}.iso_seq.fastq.gz
#        samtools fastq -@ 10 ${INDIR}${sample}.iso_seq.removed.bam | pigz -p 10 > ${sample}.iso_seq.removed.fastq.gz

        # Align reads to Cur reference genome
        minimap2 -ax splice:hq -t 20 -R '@RG\tID:${sample}\tSM:good_prim' $refcur ${sample}.iso_seq.fastq.gz \
            | samtools sort -O BAM -@ 20 \
            > ${sample}.iso_seq.bam
        samtools index ${sample}.iso_seq.bam
        bamCoverage --numberOfProcessors 20 -b ${sample}.iso_seq.bam -of bedgraph -o ${sample}.iso_seq.bedgraph
        minimap2 -ax splice:hq -t 20 -R '@RG\tID:${sample}\tSM:bad_prim' $refcur ${sample}.iso_seq.removed.fastq.gz \
            | samtools sort -O BAM -@ 20 \
            > ${sample}.iso_seq.removed.bam
        samtools index ${sample}.iso_seq.removed.bam
        bamCoverage --numberOfProcessors 20 -b ${sample}.iso_seq.removed.bam -of bedgraph -o ${sample}.iso_seq.removed.bedgraph
    "
done

# Get readcount at candidate-indel positions
canind='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/all_good_candidate_indels.qualfilter.exons.txt'
for sample in ${SAMPLES[@]}; do
    cd  $CWD; cd $sample
    bam-readcount -b 30 -q 60 -w 0 -l $canind -f $refcur ${sample}.iso_seq.filtered.bam \
| awk '{printf $1" "$2" "$3" "$4; for(i=6;i<=10;i++) {n1=split($i,a,":"); printf " "a[2]};  for(i=11;i<=NF;i++) {n1=split($i,a,":"); printf " "a[1]" "a[2]}; printf "\n"}' > ${sample}.iso_seq.b30_q60.readcount &
    bam-readcount -b 30 -q 60 -w 0 -l $canind -f $refcur ${sample}.iso_seq.removed.bam \
| awk '{printf $1" "$2" "$3" "$4; for(i=6;i<=10;i++) {n1=split($i,a,":"); printf " "a[2]};  for(i=11;i<=NF;i++) {n1=split($i,a,":"); printf " "a[1]" "a[2]}; printf "\n"}' > ${sample}.iso_seq.removed.b30_q60.readcount &
done


