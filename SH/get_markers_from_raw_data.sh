#Use raw data for the four libraries and get SNP markers from them directly.
#The markers identified from here can then be compared/called in the pooled
#for each sample.

cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/get_markers_from_raw_data/'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/bigdata/'
curidx='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta.mm2_Xsr.idx'
trimmer='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/trim_barcodes_from_scdna_reads.py'
samples=( "MUT_11_1" "MUT_15" "WT_1" "WT_19" )
bcsize=16

cd $cwd
for sample in ${samples[@]}; do
  cd $cwd
  mkdir $sample
  cd $sample
  bsub -q multicore20 -n 20 -R "span[hosts=1] rusage[mem=30000]" -M 35000 -o a.log -e b.err "
#    cat ${indir}/${sample}/L50*R1* > ${sample}_R1.fastq.gz
#    cat ${indir}/${sample}/L50*R2* > ${sample}_R2.fastq.gz
#    python $trimmer $bcsize  ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz
#    minimap2 -ax sr -t 8 $curidx ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz \
#      | samtools view -O BAM - \
#      | samtools sort - > ${sample}.sorted.bam
#      samtools index ${sample}.sorted.bam

    samtools sort -n -@ 18 -O BAM ${sample}.sorted.bam \
    | samtools fixmate -c -m -O BAM -@ 18 - - \
    | samtools sort -O BAM -@ 18 - \
    | samtools markdup -r -S -s -O BAM -@ 18 - - \
    | samtools sort -O BAM -@ 18 - > ${sample}.rmdup.sorted.bam
    samtools index ${sample}.rmdup.sorted.bam
  "
done
