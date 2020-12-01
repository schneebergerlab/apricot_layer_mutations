curidx='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta.mm2_Xsr.idx'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/rp_leaf_illumina/'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/rojo_passion/'


# align individual samples and remove duplicates
samples=( 'wt_A' 'wt_B' 'mut_C' 'mut_D' )
cd $cwd
for sample in ${samples[@]}; do
  cd $cwd
  cd $sample
  bsub -q multicore20 -n 20 -R "span[hosts=1] rusage[mem=15000]" -M 20000 -oo $sample.log -eo $sample.err "
#    skewer -r 0.1 -d 0.05 -k 8 -q 20 -l 75 -m pe -t 20 -x \
#      srv/netscratch/dep_mercier/grp_schneeberger/projects/hyper_co/data/reads/shqMergedAdapters_Primers_representative_rc.fa \
#      ${indir}/${sample}_R1.fastq.gz \
#      ${indir}/${sample}_R2.fastq.gz \
#      -z -o ${sample}_ql

    minimap2 -ax sr --eqx -t 20 $curidx \
      ${sample}_ql-trimmed-pair1.fastq.gz \
      ${sample}_ql-trimmed-pair2.fastq.gz \
    | samtools view -b - \
    | samtools sort -O BAM -o ${sample}.sorted.bam -
    samtools index ${sample}.sorted.bam

    ## Mark duplicated reads
    samtools sort -@ 10 -n -O BAM ${sample}.sorted.bam \
    | samtools fixmate -@ 10 -c -m -O BAM - - \
    | samtools sort -@ 10 -O BAM - \
    | samtools markdup -@ 10 -S -s -O BAM - - \
    | samtools sort -@ 10 -O BAM - > ${sample}.DUPmarked.bam
    samtools index ${sample}.DUPmarked.bam

    samtools view -G 1024 -O BAM ${sample}.DUPmarked.bam > ${sample}.DUPmarked.deduped.bam
    samtools index ${sample}.DUPmarked.deduped.bam
  "
done

# Merge samples to get representative tree sequence
cd $cwd
samtools merge -O BAM -@ 60 merged_samples/merged.bam */*.DUPmarked.deduped.bam
cd merged_samples
samtools index merged.bam