cur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta.mm2_Xsr.idx'
cur_contig_list='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta.contig_list'
samples=( 'currot' 'orangered' )
#samples=( 'currot' )
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/'
ploidy=/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Plum/results/marker_generated/ploidy

mq=10
for sample in ${samples[@]}; do
  cd ${cwd}/$sample
  bsub -q multicore40 -n 40  -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo align_${sample}.log -eo align_${sample}.err "
    minimap2 -ax sr -t 40 $cur *R1*.fastq.gz *R2*.fastq.gz \
    | samtools view -O BAM - \
    | samtools sort - > ${sample}_cur_v1.1.sorted.bam
    samtools index ${sample}_cur_v1.1.sorted.bam

    ## Remove duplicates
    samtools sort -n -@ 38 -O BAM ${sample}_cur_v1.1.sorted.bam \
    | samtools fixmate -c -m -O BAM -@ 38 - - \
    | samtools sort -O BAM -@ 38 - \
    | samtools markdup -r -S -s -O BAM -@ 38 - - \
    | samtools sort -O BAM -@ 38 - > ${sample}_cur_v1.1.rmdup.sorted.bam
    samtools index ${sample}_cur_v1.1.rmdup.sorted.bam

    ## Get readcount
    bam-readcount -b 30 -q 20 -w 0 -f $refcur ${sample}_cur_v1.1.rmdup.sorted.bam| awk '{if($4>3) {n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":"); n6=split($11,f, ":"); n7=split($12,g, ":"); n8=split($13,h, ":");  print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2], f[1], f[2], g[1], g[2], h[1], h[2]}}' > bam_read_counts_b30_q20.txt

    ## Get positions with non-reference alleles
    python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py bam_read_counts_b30_q20.txt

    samtools view -q 10 -@ 38 -b ${sample}_cur_v1.1.rmdup.sorted.bam > ${sample}_cur_v1.1_q10.rmdup.sorted.bam
    samtools index ${sample}_cur_v1.1_q10.rmdup.sorted.bam

    cat $cur_contig_list \
    | sed 's/\\n/ /g' \
    | xargs -d '\n' -n 1 -P 40 -I {} /bin/bash -c ' \
        bcftools mpileup -q \$2 -d 800 --threads 1 -f \$3 -r \$1 \${4}.rmdup.sorted.bam  > \${4}_mq\${2}_\${1}.pileup \
      ' -- {} $mq cur ${sample}_cur_v1.1_q10

    bcftools concat --threads 38 ${sample}_cur_v1.1_q10_mq${mq}_*pileup \
    | bcftools call --threads 38 -m -v -Ov --ploidy-file ${ploidy} \
    > ${sample}_cur_v1.1_q10_mq${mq}.vcf

    vcftools --vcf ${sample}_cur_v1.1_q10_mq${mq}.vcf --freq --out ${sample}_cur_v1.1_q10_mq${mq}
    vcftools --vcf ${sample}_cur_v1.1_q10_mq${mq}.vcf --remove-indels --recode --recode-INFO-all --out ${sample}_cur_v1.1_q10_onlySNP_mq${mq}.vcf
  "
done


