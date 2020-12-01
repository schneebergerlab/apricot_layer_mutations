#Use raw data for the four libraries and get SNP markers from them directly.
#The markers identified from here can then be compared/called in the pooled
#for each sample.

cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/get_markers_from_raw_data/'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/bigdata/'
curidx='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta.mm2_Xsr.idx'
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
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

#
#    ## Mark duplicated reads
#    samtools sort -n -O BAM ${sample}.sorted.bam \
#    | samtools fixmate -c -m -O BAM - - \
#    | samtools sort -O BAM - \
#    | samtools markdup -S -s -O BAM - - \
#    | samtools sort -O BAM - > ${sample}.DUPmarked.bam
#    samtools index ${sample}.DUPmarked.bam

#    # Get non-dup reads in fastq.gz from the markdup bam files
    python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/remove_duplicates_from_marked_bam_file.py ${sample}.DUPmarked.bam -p ${sample}.dedup

    gzip ${sample}.dedup_R1.fastq
    gzip ${sample}.dedup_R2.fastq

     minimap2 -ax sr -t 17 $curidx ${sample}.dedup_R1.fastq.gz ${sample}.dedup_R2.fastq.gz \
      | samtools view -O BAM - \
      | samtools sort - > ${sample}.dedup.sorted.bam
      samtools index ${sample}.dedup.sorted.bam

    bam-readcount -b 30 -q 20 -w 0 -f $refcur ${sample}.dedup.sorted.bam | awk '{if(\$4>3) {n1=split(\$6,a,\":\"); n2=split(\$7,b,\":\"); n3=split(\$8,c, \":\"); n4=split(\$9,d,\":\"); n5=split(\$10,e,\":\"); n6=split(\$11,f, \":\"); n7=split(\$12,g, \":\"); n8=split(\$13,h, \":\");  print \$1, \$2, \$3, \$4, a[2], b[2], c[2], d[2], e[2], f[1], f[2], g[1], g[2], h[1], h[2]}}' > bam_read_counts_b30_q20.txt

  # GET POSITIONS WITH AT LEAST THREE NON-REFERENCE BASES
  python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py bam_read_counts_b30_q20.txt
  "
done
