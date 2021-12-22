# Link data to working folder
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_isoseq/
cd $CWD
ln -s /srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4876/4876_A_run545_CCS.bam WT_1.iso_seq.ccs.bam
ln -s /srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4876/4876_B_run545_CCS.bam WT_19.iso_seq.ccs.bam
ln -s /srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4876/4876_C_run545_CCS.bam MUT_11_1.iso_seq.ccs.bam
ln -s /srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4954/4954_A_run545_CCS.bam MUT_15.iso_seq.ccs.bam

################################################################################
############ Clusgter Iso-seq reads into individual barcodes ###################
################################################################################

## Select reads having good primer orientation
### Lima from syri3.8 environment
SAMPLES=("MUT_11_1" "MUT_15" "WT_1" "WT_19")
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_isoseq/
for sample in ${SAMPLES[@]}; do
    lima --isoseq --dump-clips --dump-removed ${sample}.iso_seq.ccs.bam 10x_primers.fasta ${sample}.iso_seq.bam &
done

## Get BC and UMI tags for the reads
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    isoseq3 tag ${s}.iso_seq.5p--3p.bam ${s}.iso_seq.fl.bam --design T-12U-16B &
done
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    isoseq3 refine --require-polya ${s}.iso_seq.fl.bam 10x_primers.fasta ${s}.iso_seq.flnc.bam &
done


### Select reads with correct primer orientation from the raw CCS bam
#for s in WT_1 WT_19 MUT_11_1 MUT_15; do
#    {
#    samtools view -@ 12 ${s}.iso_seq.5p--3p.bam | cut -f1 > ${s}_read_names.txt
#    samtools view -@ 12 -H ${s}.iso_seq.ccs.bam > ${s}.iso_seq.ccs.good_primer.sam
#    samtools view -@ 12 ${s}.iso_seq.ccs.bam \
#    | grep -Ff ${s}_read_names.txt \
#    >> ${s}.iso_seq.ccs.good_primer.sam
#    samtools view -@ 12 -b ${s}.iso_seq.ccs.good_primer.sam > ${s}.iso_seq.ccs.good_primer.bam
#    rm ${s}.iso_seq.ccs.good_primer.sam
#    } &
#done
#
### Rerun Lima but this time append the BCs (from scRNA) to the 10X_primers.fasta
#### Append the BC to primers
#BCDIR=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/barcodes/
#for s in wt1 wt19 mut11 mut15; do
#    cd $CWD
#    rm ${s}_bcs.fasta
#    while read r; do
#        bc=$(echo $r | cut -d' ' -f 2)
#        echo -e ">${bc}_5p\nCTACACGACGCTCTTCCGATCT${bc}\n" >> ${s}_bcs.fasta
#    done < $BCDIR/${s}_bcs_clstr_id.txt
#    echo -e ">${end}_3p\nAAGCAGTGGTATCAACGCAGAGTACATGGG\n" >> ${s}_bcs.fasta
#done
#mv wt1_bcs.fasta WT_1_bcs.fasta
#mv wt19_bcs.fasta WT_19_bcs.fasta
#mv mut11_bcs.fasta MUT_11_1_bcs.fasta
#mv mut15_bcs.fasta MUT_15_bcs.fasta
#for s in WT_1 WT_19 MUT_11_1 MUT_15; do
#    cd $CWD
#    mkdir ${s}_no3primer; cd ${s}_no3primer
#    {
#    lima --isoseq -d \
#    --dump-clips \
#    --dump-removed \
#    -j 12 \
#    ../${s}.iso_seq.ccs.good_primer.bam \
#    ../${s}_bcs.no3primer.fasta \
#    ${s}.iso_seq.bam
#    } &
#done
#
#### Test without using the 3' primer
#for s in wt1 wt19 mut11 mut15; do
#    cd $CWD
#    rm ${s}_bcs.no3primer.fasta
#    while read r; do
#        bc=$(echo $r | cut -d' ' -f 2)
#        echo -e ">${bc}_5p\n${bc}\n" >> ${s}_bcs.no3primer.fasta
#    done < $BCDIR/${s}_bcs_clstr_id.txt
#    echo -e ">${end}_3p\nAAGCAGTGGTATCAACGCAGAGTACATGGG\n" >> ${s}_bcs.no3primer.fasta
#done
#mv wt1_bcs.no3primer.fasta WT_1_bcs.no3primer.fasta
#mv wt19_bcs.no3primer.fasta WT_19_bcs.no3primer.fasta
#mv mut11_bcs.no3primer.fasta MUT_11_1_bcs.no3primer.fasta
#mv mut15_bcs.no3primer.fasta MUT_15_bcs.no3primer.fasta
#for s in WT_1 WT_19 MUT_11_1 MUT_15; do
#    cd $CWD
#    mkdir ${s}_no3primer; cd ${s}_no3primer
#    {
#    lima --isoseq -d \
#    --dump-clips \
#    --dump-removed \
#    -j 12 \
#    ../${s}.iso_seq.ccs.good_primer.bam \
#    ../${s}_bcs.no3primer.fasta \
#    ${s}.iso_seq.bam
#    } &
#done
##### Without the 3' primer fewer reads were classified into BCs

## Get reads scRNA BCs
INDIR='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_isoseq/'
CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/'
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $CWD
    mkdir $s; cd $s
    samtools sort -t XC -O BAM ${INDIR}/${s}.iso_seq.flnc.bam > ${s}.XC_sorted.bam &
    samtools index ${s}.XC_sorted.bam
    samtools view $s.XC_sorted.bam  | grep -o -P 'XC:Z:[^\s]*' | uniq -c > cnt.txt &
done
/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/iso_seq_analysis.py
cd $CWD
pigz -p 10 */*_reads.fa

## Map reads
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/'
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $CWD; cd $s
    clst=$(ls clstrs*reads.fa.gz)
    for c in ${clst[@]}; do
        cls=$(echo $c | sed 's/_reads.fa.gz//g')
        bsub -q multicore20 -n 10 -R "span[hosts=1] rusage[mem=20000]" -M 25000 -oo ${cls}.log -eo ${cls}.err "
            minimap2 -ax splice:hq -t 10 -R '@RG\tID:${cls}\tSM:1' --secondary=no -uf -C5 -O6,24 -B4 -Y $refcur $c \
            | samtools view -F 2048 -h - \
            | samtools sort -O BAM -@ 10 - \
            > ${cls}.iso_seq.bam
            samtools index -@ 10 ${cls}.iso_seq.bam
    "
    done
done

## Get readcount at candidate-indel positions
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
CWD='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/'
muts='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mutations.bed'
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $CWD; cd $s
    clstrs=$(ls clstrs*bam)
    for cls in ${clstrs[@]}; do
        c=$(echo $cls | sed 's/bam/rc.txt/g')

        bam-readcount -b 0 -q 0 -w 0 -l $muts -f $refcur $cls \
        | awk '{printf $1" "$2" "$3" "$4; for(i=6;i<=10;i++) {n1=split($i,a,":"); printf " "a[2]};  for(i=11;i<=NF;i++) {n1=split($i,a,":"); printf " "a[1]" "a[2]}; printf "\n"}' > $c &
#        hometools pbamrc \
#        -l $muts \
#        -f $refcur \
#        -w 0 \
#        -n 1 \
#        $cls \
#        $c \
#        &
    done
done



for sample in ${SAMPLES[@]}; do
    cd  $CWD; cd $sample
    bam-readcount -b 30 -q 60 -w 0 -l $canind -f $refcur ${sample}.iso_seq.filtered.bam \
| awk '{printf $1" "$2" "$3" "$4; for(i=6;i<=10;i++) {n1=split($i,a,":"); printf " "a[2]};  for(i=11;i<=NF;i++) {n1=split($i,a,":"); printf " "a[1]" "a[2]}; printf "\n"}' > ${sample}.iso_seq.b30_q60.readcount &
    bam-readcount -b 30 -q 60 -w 0 -l $canind -f $refcur ${sample}.iso_seq.removed.bam \
| awk '{printf $1" "$2" "$3" "$4; for(i=6;i<=10;i++) {n1=split($i,a,":"); printf " "a[2]};  for(i=11;i<=NF;i++) {n1=split($i,a,":"); printf " "a[1]" "a[2]}; printf "\n"}' > ${sample}.iso_seq.removed.b30_q60.readcount &
done


