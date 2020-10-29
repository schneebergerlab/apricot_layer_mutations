cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/ryu_et_al_plant_physio_2019/'

athal_ref_idx='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/athal/TAIR10_chr_all.fa.mm2_Xsr.idx'
athal_ref='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/athal/TAIR10_chr_all.fa'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/ryu_et_al_plant_physio_2019/reads/'
gene_bed='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/athal/geneCoords.bed'
cd $cwd
#fqs=( 'wt1.fastq.gz'	'wt2.fastq.gz'  'wt3.fastq.gz'  'gl2.fastq.gz' 'rhd6.fastq.gz'	)
#############################################################
## Step 1: Index Reference Genome
#############################################################
#STAR  --runThreadN 40 \
#      --runMode genomeGenerate \
#      --genomeDir /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/athal/star_index \
#      --genomeFastaFiles /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/athal/TAIR10_chr_all.fa \
#      --sjdbGTFfile /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/athal/TAIR10.gtf \
#      --sjdbOverhang 149         # We do not need this. Need to rerun the indexing and alignment

#############################################################
## Step 2: Align reads to genome index
#############################################################X
ref_idx='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/athal/star_index/'
fqs=( 'wt1'	'wt2'  'wt3'  'gl2' 'rhd6'	)
for fq in ${fqs[@]}; do
  bsub -q multicore20 -n 20 -R "span[hosts=1] rusage[mem=20000]" -M 20000 -oo ${fq}.log -eo ${fq}.err "
    STAR --runThreadN 20 \
        --genomeDir $ref_idx \
        --readFilesIn ${indir}${fq}.fastq.gz \
        --readFilesCommand zcat \
        --outStd BAM_SortedByCoordinate \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMunmapped Within \
        --outSAMprimaryFlag AllBestScore \
        --outFileNamePrefix ${cwd}${fq} \
        --quantMode GeneCounts \
        --twopassMode Basic \
        --alignIntronMax 20000 \
        --outFilterMismatchNmax 5 > ${fq}.sorted.bam
    samtools index ${fq}.sorted.bam
    samtools depth -b $gene_bed -d 0 --reference $athal_ref -a ${fq}.sorted.bam > ${fq}.depth
  "
done

