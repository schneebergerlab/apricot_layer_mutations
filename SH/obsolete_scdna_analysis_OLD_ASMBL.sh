## Step 1: Index current Currot and OrangeRed assemblies to be used as reference genomes for barcode correction
<<Useful
# Find Repeatitive sequences
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/repeat_region/
BuildDatabase -name currot.v1.1 ../currot.v1.1.fasta
RepeatModeler -database currot.v1.1 -pa 60 -LTRStruct

# Find repeatitive regions using RepeatMasker
cd ..
RepeatMasker \
  -pa 60 \
  -s \
  -norna \
  -lib ./repeat_region/currot.v1.1-families.fa \
  -frag 100000 \
  -dir ./repeat_region/ \
  -poly \
  -u currot.v1.1.fasta

tail +4  currot.v1.1.fasta.out |  awk '{print $5"\t"$6-1"\t"$7}'  > currot.v1.1.fasta.out.bed

## Add buffer of 20 BP around indel positions
cat currot.v1.1.fasta.out.bed \
| awk '{printf $1"\t"} {if($2<20) {printf 0"\t"} else {printf $2-20"\t"}} {print $3+20}' \
> currot.v1.1.fasta.out.padded_20.bed
Useful


<<Useful
indir=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/

cd $indir
cellranger-dna mkref currot.v1.1.fasta currot_contig_defs.json
Useful

## Step 2: Use the indexed genomes and run the cellranger-dna CNV identification pipeline. One of the internediate steps for this pipeline is to barcode correction. We can use the final bam file to get the corrected barcodes.
## Currot genome would be used as the reference genome for all analysis (preliminary testing showed that there are not major differences between currot vs orangered)

run_barcode_correction () {
    bsub -q ioheavy  -n 40 -R "span[hosts=1] rusage[mem=100000]" -M 100000 -oo $1.log -eo $1.err "
    export _JAVA_OPTIONS=-Xgcthreads4
        cellranger-dna cnv --id=$1 --reference=$2 --fastq=$3 --sample=$4 --localcores=40 --localmem=98
        "
}
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/refdata-currot.v1.1/'
indir="/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/bigdata/"
python="/netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/syri3.7/bin/python"
filter="/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/filter_short_reads.py"

samples=( "MUT_11_1" "MUT_15" "WT_1" "WT_19" )

<<Useful
##  All reads smaller than 50bps were filtered out as cellranger-dna cannot work with smaller reads and there were naked barcode sequenced as reads. (Question: why did naked barcode ended up getting sequenced?)
for sample in ${samples[@]}; do
  cd ${indir}/${sample}/
  R1s=($(ls ${sample}*R1*))
  R2s=($(ls ${sample}*R2*))
  cnt=$(echo ${R2s[@]} | wc -w)
  for ((i=0;i<cnt;i++)); do
    bsub -q normal -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo ${i}.log -eo ${i}.err "
      $python $filter ${R1s[${i}]} ${R2s[$i]} 50
    "
  done
done
Useful

## Cellranger pipeline on all samples
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/bigdata/'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/cellranger_out'
<<Useful
for sample in ${samples[@]}; do
  cd $cwd
  R2s=($(ls ${indir}/${sample}/renamed_L50_${sample}*R2* | xargs -n 1 basename))
  pres=''
  for R2 in ${R2s[@]}; do
    pres=$(echo ${pres}${R2%_S*},)
  done
  pres=$(echo ${pres} | sed 's/,$//g')
  run_barcode_correction $sample $refcur ${indir}$sample/ $pres
done
Useful

## Step 3: Separate reads from the corrected barcode bam file to separate folders corresponding to each barcode, get read count, and align them to the reference
T10xbam2fq=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/binaries/T10xbam2fq
asCellseparator=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/binaries/asCellseparator
curidx=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta.mm2_Xsr.idx

MIN_RC=10000

get_cell_reads () {
    bsub -q multicore20 -n 20 -J ${2}_getbc -R "span[hosts=1] rusage[mem=20000]" -M 20000 -m "hpc001.cluster.loc hpc002.cluster.loc" -oo get_cell_reads.log -eo get_cell_reads.err "
            samtools sort -@ 19 -n $1 -o ${2}_RNsorted_bam.bam
            samtools view ${2}_RNsorted_bam.bam | $T10xbam2fq - ${2}
            mkdir barcodes
            $asCellseparator 16 ${2}_fqfrom10xBam_bxCorrected_R1.fq.gz ${2}_fqfrom10xBam_bxCorrected_R2.fq.gz $MIN_RC 10000000 barcodes/
            "
}

merge_fastqs () {
#    bsub -q normal -w "done(${1}_getbc)" -J ${1}_merge -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo merge_fastqs.log -eo merge_fastqs.err "
    bsub -q normal -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo merge_fastqs.log -eo merge_fastqs.err "
        readlink -f barcodes/* > barcodes_list
        while read r; do
            cd \$r
            bc=\$(basename \$r)
            cat part*R1*.fq.gz > \${bc}_R1.fastq.gz
            cat part*R2*.fq.gz > \${bc}_R2.fastq.gz
        done < barcodes_list

        cd ${2}/${1}
        rm barcodes_read_count
        while read r; do
            bc=\$(basename \$r)
            cnt=\$(zgrep -c '@' \$r/\${bc}_R1.fastq.gz)
            echo \$bc \$cnt >> barcodes_read_count
        done < barcodes_list
        "
}

align_cells () {
#    bsub -q multicore40 -n 20 -w "done(${1}_merge)" -J ${1}_align -R "span[hosts=1] rusage[mem=25000]" -M 30000 -oo align_cells.log -eo align_cells.err "
    bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=60000]" -M 75000 -oo align_cells.log -eo align_cells.err "
    xargs -a barcodes_list -n 1 -P 40 -I {} /bin/bash -c '
        cd {}
        bc=\$(basename {})
#        minimap2 -ax sr -t 1 \$2 \${bc}_R1.fastq.gz \${bc}_R2.fastq.gz \
#        | samtools view -O BAM - \
#        | samtools sort - > \${bc}.sorted.bam
#        samtools index \${bc}.sorted.bam

        ## Mark duplicated reads
#        samtools sort -n -O BAM \${bc}.sorted.bam \
#        | samtools fixmate -c -m -O BAM - - \
#        | samtools sort -O BAM - \
#        | samtools markdup -S -s -O BAM - - \
#        | samtools sort -O BAM - > \${bc}.DUPmarked.bam
#        samtools index \${bc}.DUPmarked.bam

        # Get non-dup reads in fastq.gz from the markdup bam files
        python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/remove_duplicates_from_marked_bam_file.py \${bc}.DUPmarked.bam -p \${bc}_dedup

        gzip -f \${bc}_dedup_R1.fastq
        gzip -f \${bc}_dedup_R2.fastq

        # Get non-dup bam file from the markdup bam file
#        samtools view -G 1024 -O BAM \${bc}.DUPmarked.bam > \${bc}.DUPmarked.deduped.bam
#        samtools index \${bc}.DUPmarked.deduped.bam
#
#        # Get BAM coverage
#        bamCoverage --numberOfProcessors 1 -b \${bc}.DUPmarked.deduped.bam -of bedgraph -o \${bc}.bedgraph
' -- {} $2
    "
}

samples=( "MUT_11_1" "MUT_15" "WT_1" "WT_19" )
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/cellranger_out/'

<<Useful

for sample in ${samples[@]}; do
    cd $cwd/
    mkdir $sample
    cd $sample
#    get_cell_reads ${indir}/${sample}/outs/possorted_bam.bam $sample
#    merge_fastqs $sample $cwd
    align_cells $sample $curidx
done

Useful

####################################################################
############ Step  3: Aligned pooled data reads and get read-counts at variant positions
####################################################################
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/'
curidx='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta.mm2_Xsr.idx'
curidxbt2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1'
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
samples=( "MUT_11_1" "MUT_15" "WT_1" "WT_19" )

cd $cwd

## Step 3a: Merge individual barcodes fastq to one sample fastq
for sample in ${samples[@]}; do
  cd $cwd
  mkdir $sample
  cd $sample
  rm ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz
  while read r; do
    bc=$(basename $r)
#    cat ${indir}${sample}/barcodes/${bc}/${bc}_dedup_R1.fastq | gzip -c >> ${sample}_R1.fastq.gz
#    cat ${indir}${sample}/barcodes/${bc}/${bc}_dedup_R2.fastq | gzip -c >> ${sample}_R2.fastq.gz
    cat ${indir}${sample}/barcodes/${bc}/${bc}_dedup_R1.fastq.gz  >> ${sample}_R1.fastq.gz
    cat ${indir}${sample}/barcodes/${bc}/${bc}_dedup_R2.fastq.gz  >> ${sample}_R2.fastq.gz
  done < ${indir}${sample}/barcodes_list
done

## Step 3b: Align merged reads against the currot reference genome and get potential variant positions using minimap2
for sample in ${samples[@]}; do
  cd $cwd
  mkdir $sample
  cd $sample
  bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo align_reads_to_cur.log -eo align_reads_to_cur.err "
    minimap2 -ax sr -t 30 $curidx ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz \
    | samtools view -O BAM - \
    | samtools sort - > ${sample}.sorted.bam
    samtools index ${sample}.sorted.bam

#     Duplicates were removed from the individual cells.
#     Therefore, there is no need to remove duplicates from the merged data.
#     All 'duplicated' reads in the merged data would corresponds to actual resequencing of the region

    bam-readcount -b 30 -q 20 -w 0 -f $refcur ${sample}.sorted.bam | awk '{if(\$4>3) {n1=split(\$6,a,\":\"); n2=split(\$7,b,\":\"); n3=split(\$8,c, \":\"); n4=split(\$9,d,\":\"); n5=split(\$10,e,\":\"); n6=split(\$11,f, \":\"); n7=split(\$12,g, \":\"); n8=split(\$13,h, \":\");  print \$1, \$2, \$3, \$4, a[2], b[2], c[2], d[2], e[2], f[1], f[2], g[1], g[2], h[1], h[2]}}' > bam_read_counts_b30_q20.txt

  # GET POSITIONS WITH AT LEAST THREE NON-REFERENCE BASES
  python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py bam_read_counts_b30_q20.txt
  "
done

## Step 3c: Align merged reads against the currot reference genome and get potential variant positions using bowtie2
for sample in ${samples[@]}; do
  cd $cwd
  mkdir $sample
  cd $sample
  bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo align_reads_to_cur.log -eo align_reads_to_cur.err "
    bowtie2 --end-to-end \
            --very-sensitive \
            --threads 40 \
            -x $curidxbt2 \
            -1 ${sample}_R1.fastq.gz \
            -2 ${sample}_R2.fastq.gz \
            -S ${sample}.bt2.sam

    samtools view -O BAM ${sample}.bt2.sam \
    | samtools sort - > ${sample}.sorted.bt2.bam
    samtools index ${sample}.sorted.bt2.bam

    bam-readcount -b 30 -q 20 -w 0 -f $refcur ${sample}.sorted.bt2.bam | awk '{if(\$4>3) {n1=split(\$6,a,\":\"); n2=split(\$7,b,\":\"); n3=split(\$8,c, \":\"); n4=split(\$9,d,\":\"); n5=split(\$10,e,\":\"); n6=split(\$11,f, \":\"); n7=split(\$12,g, \":\"); n8=split(\$13,h, \":\");  print \$1, \$2, \$3, \$4, a[2], b[2], c[2], d[2], e[2], f[1], f[2], g[1], g[2], h[1], h[2]}}' > bam_read_counts_b30_q20.bt2.txt

  # GET POSITIONS WITH AT LEAST THREE NON-REFERENCE BASES
  python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py bam_read_counts_b30_q20.bt2.txt
  "
done

####################################################################
############ Step 4: Get mutated positions
####################################################################
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
test='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_1/head.bam'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
samples=( "MUT_11_1" "MUT_15" "WT_1" "WT_19" )

## Step 4a: Get candidate variant positions using heuristic cut-offs
cd $cwd
## Use python3.7 environment
python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_candidates_variant_positions.py \
 WT_1/filtered_low_ref_al_bam_read_counts_b30_q20.txt \
 WT_19/filtered_low_ref_al_bam_read_counts_b30_q20.txt \
 MUT_11_1/filtered_low_ref_al_bam_read_counts_b30_q20.txt \
 MUT_15/filtered_low_ref_al_bam_read_counts_b30_q20.txt \
 -s WT_1 WT_19 MUT_11_1 MUT_15 &

python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_candidates_variant_positions.py \
 WT_1/filtered_low_ref_al_bam_read_counts_b30_q20.bt2.txt \
 WT_19/filtered_low_ref_al_bam_read_counts_b30_q20.bt2.txt \
 MUT_11_1/filtered_low_ref_al_bam_read_counts_b30_q20.bt2.txt \
 MUT_15/filtered_low_ref_al_bam_read_counts_b30_q20.bt2.txt \
 -s WT_1_bt2 WT_19_bt2 MUT_11_1_bt2 MUT_15_bt2 &

for sample in ${samples[@]}; do
  cd $cwd
  cd $sample
  bedtools intersect -a ${sample}_only_SNPs_candidate.sorted.bed -b ${sample}_bt2_only_SNPs_candidate.sorted.bed > ${sample}_only_SNPs_candidate.sorted.common.bed
  cat ${sample}_only_SNPs_candidate.sorted.common.bed | awk '{print $1"\t"$2+1"\t"$3"\t"$4"\t"$5}' > ${sample}_only_SNPs_candidate.sorted.common.regions
done

## Step 4b: Select candidates supported by multiple cells
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
samples=( "MUT_11_1" "MUT_15" "WT_1" "WT_19" )
for sample in ${samples[@]}; do
  cd $cwd
  cd $sample
  barcodes_list=${indir}${sample}/barcodes_list
  bedfile=$(rf ${sample}_only_SNPs_candidate.sorted.common.regions)
  inpath=${indir}${sample}/barcodes/
  mkdir cells_readcount
  cd cells_readcount
  while read r; do
    bc=$(basename $r)
    mkdir $bc
    cd $bc
    bsub -q short -R "span[hosts=1] rusage[mem=2000]" -M 2500 -oo rc_at_candidate.log -eo rc_at_candidate.err "
      bam-readcount -b 30 -q 20 -w 0 -l $bedfile -f $refcur ${inpath}/${bc}/${bc}.DUPmarked.deduped.bam \
      | awk '{n1=split(\$6,a,\":\"); n2=split(\$7,b,\":\"); n3=split(\$8,c, \":\"); n4=split(\$9,d,\":\"); n5=split(\$10,e,\":\"); print \$1, \$2, \$3, \$4, a[2], b[2], c[2], d[2], e[2]}' > ${bc}_candidate_readcount.txt
    "
    cd ..
  done < $barcodes_list
done

## Use python3.7 environment
for sample in ${samples[@]}; do
  cd $cwd
  cd $sample
  echo $sample
  python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/select_candidate_supported_by_cells.py \
  ${sample}_only_SNPs_candidate.txt \
  ${sample}_bt2_only_SNPs_candidate.txt \
  ${sample}_only_SNPs_candidate.sorted.common.regions \
  cells_readcount/
done


## Step 4c: Select candidates that are in highly unqiue (low mappability) regions
# Use genmap to get regions which have low mappability (high uniqueness) in the genomes
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose
genmap index -F currot.v1.1.fasta -I currot.v1.1.fasta.genmap.idx
genmap map -I currot.v1.1.fasta.genmap.idx/ -O currot.v1.1.fasta.genmap.E0_K31.map -E 0 -K 31 -w -bg -T 60
cat currot.v1.1.fasta.genmap.E0_K31.map.bedgraph | awk '{if($4==1) print $0}' > currot.v1.1.fasta.genmap.E0_K31.map.unique.bedgraph

genmap map -I currot.v1.1.fasta.genmap.idx/ -O currot.v1.1.fasta.genmap.E0_K51.map -E 0 -K 51 -w -bg -T 60
cat currot.v1.1.fasta.genmap.E0_K51.map.bedgraph | awk '{if($4==1) print $0}' > currot.v1.1.fasta.genmap.E0_K51.map.unique.bedgraph


cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
inbed='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta.genmap.E0_K51.map.unique.bedgraph'
samples=( "MUT_11_1" "MUT_15" "WT_1" "WT_19" )
cd $cwd
## Use python3.7 environment
for sample in ${samples[@]}; do
  cd $cwd
  cd $sample
  echo $sample

  # Select candidates that are in uniquely mapping regions
  cat multi_cell_${sample}_only_SNPs_candidate.sorted.common.regions \
  | awk '{if($2>=26 && $3>=26){print $1"\t"$2-26"\t"$3-25"\t"$4"\t"$5"\t"$6"\t"$7}}' \
  > multi_cell_${sample}_only_SNPs_candidate.sorted.common.regions.adjusted.bed
  bedtools intersect \
  -a multi_cell_${sample}_only_SNPs_candidate.sorted.common.regions.adjusted.bed \
  -b $inbed \
  | awk '{print $1"\t"$2+26"\t"$3+25"\t"$4"\t"$5"\t"$6"\t"$7}' \
  > multi_cell_${sample}_only_SNPs_candidate.sorted.common.unique_genome.regions
#  rm multi_cell_${sample}_only_SNPs_candidate.sorted.common.regions.adjusted.bed

  # Remove candidates near repeatitive regions
  cat multi_cell_${sample}_only_SNPs_candidate.sorted.common.unique_genome.regions \
  | awk '{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
  > multi_cell_${sample}_only_SNPs_candidate.sorted.common.unique_genome.bed
  bedtools intersect \
    -v \
    -a multi_cell_${sample}_only_SNPs_candidate.sorted.common.unique_genome.bed \
    -b /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/repeat_region/currot.v1.1.fasta.out.padded_20.bed \
  | awk '{print $1"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
  > multi_cell_${sample}_only_SNPs_candidate.sorted.common.unique_genome.non_repeat.regions


  # Remove candidates near indels
  cat multi_cell_${sample}_only_SNPs_candidate.sorted.common.unique_genome.non_repeat.regions \
  | awk '{print $1"\t"$2-1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
  > multi_cell_${sample}_only_SNPs_candidate.sorted.common.unique_genome.non_repeat.bed
  bedtools intersect \
    -v \
    -a multi_cell_${sample}_only_SNPs_candidate.sorted.common.unique_genome.non_repeat.bed \
    -b indels/deletions.sorted.merged.padded_20.bed \
    > no_del.bed
  bedtools intersect \
  -v \
  -a no_del.bed \
  -b indels/insertions.sorted.merged.padded_20.bed \
  | awk '{print $1"\t"$2+1"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
  > multi_cell_${sample}_only_SNPs_candidate.sorted.common.unique_genome.non_repeat.non_indel.regions
#  rm multi_cell_${sample}_only_SNPs_candidate.sorted.common.unique_genome.bed

done

# Step 4d: Perform read-level filtering


####################################################################
############ Step 5: Tool based sSNV variant calling
####################################################################

## STEP 5-Pre-requisite : Call indels using Manta, samtools, and GATK HaplotypeCaller

cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
samples=( "MUT_11_1" "MUT_15" "WT_1" "WT_19" )

# Manta
manta_config='/srv/netscratch/dep_mercier/grp_schneeberger/software/manta-1.6.0.centos6_x86_64/bin/configManta.py'
python='/netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/py2.6mg/bin/python'

cd $cwd
for sample in ${samples[@]}; do
  cd $cwd/$sample
  mkdir indels
  cd indels
  mkdir manta
  cd manta
  echo "IN $sample"
  for sample2 in ${samples[@]}; do
    if [[ $sample == $sample2 ]]; then
      continue
    else
      ## Submit MM2 based jobs
      bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${sample}_${sample2}_mm2.log -eo ${sample}_${sample2}_mm2.err "
         $python \
          $manta_config \
            --tumourBam ${cwd}/${sample}/${sample}.sorted.RG.bam \
            --normalBam ${cwd}/${sample2}/${sample2}.sorted.RG.bam \
            --referenceFasta $refcur \
            --runDir ${sample}_vs_${sample2}_mm2 \
            --scanSizeMb 2
        cd ${sample}_vs_${sample2}_mm2
        $python \
          ./runWorkflow.py \
            -m local \
            -j 40 \
            -g 50
        gunzip -c results/variants/candidateSmallIndels.vcf.gz > candidateSmallIndels.vcf
        vcf2bed --do-not-sort --insertions <candidateSmallIndels.vcf> candidateSmallIndels.insertions.bed
        vcf2bed --do-not-sort --deletions <candidateSmallIndels.vcf> candidateSmallIndels.deletions.bed

     "

      ## Submit BT2 based jobs
      bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${sample}_${sample2}_bt2.log -eo ${sample}_${sample2}_bt2.err "
        $python \
          $manta_config \
            --tumourBam ${cwd}/${sample}/${sample}.sorted.RG.bt2.bam \
            --normalBam ${cwd}/${sample2}/${sample2}.sorted.RG.bt2.bam \
            --referenceFasta $refcur \
            --runDir ${sample}_vs_${sample2}_bt2 \
            --scanSizeMb 2
        cd ${sample}_vs_${sample2}_bt2
        $python \
          ./runWorkflow.py \
            -m local \
            -j 40 \
            -g 50
        gunzip -c results/variants/candidateSmallIndels.vcf.gz > candidateSmallIndels.vcf
        vcf2bed --do-not-sort --insertions <candidateSmallIndels.vcf> candidateSmallIndels.insertions.bed
        vcf2bed --do-not-sort --deletions <candidateSmallIndels.vcf> candidateSmallIndels.deletions.bed
      "
    fi
  done
done

# Samtools
cur_contig_list='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta.contig_list'
ploidy='/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Plum/results/marker_generated/ploidy'
cd $cwd
for sample in ${samples[@]}; do
  cd $cwd/$sample
  cd indels
  mkdir samtools
  cd samtools

  # Run job for MM2 alignments
  mkdir mm2
  cd mm2
  bsub -q normal -n 1  -R "span[hosts=1] rusage[mem=5000]" -M 6000 -o samtools_mm2.log -e samtools_mm2.err "
    cat $cur_contig_list \
    | sed 's/\\n/ /g' \
    | xargs -d '\n' -n 1 -P 40 -I {} /bin/bash -c ' \
        bcftools mpileup -q \$2 -d 800 --threads 1 -f \$3 -r \$1 \${4} > \${1}_mq\${2}.pileup \
      ' -- {} 1 $refcur ${cwd}/${sample}/${sample}.sorted.RG.bam

    bcftools concat --threads 38 *_mq1.pileup \
    | bcftools call --threads 38 -m -v -Ov --ploidy-file ${ploidy} \
    > ${sample}_mq1.vcf

    vcftools --vcf ${sample}_mq1.vcf --keep-only-indels --recode --recode-INFO-all --out ${sample}_mq1_onlyindels.vcf

    vcf2bed --do-not-sort --insertions <${sample}_mq1_onlyindels.vcf.recode.vcf> ${sample}_mq1_onlyindels.vcf.recode.insertions.bed
    vcf2bed --do-not-sort --deletions <${sample}_mq1_onlyindels.vcf.recode.vcf> ${sample}_mq1_onlyindels.vcf.recode.deletions.bed

    "

  # Run job for bt2 alignments
  cd ..
  mkdir bt2
  cd bt2
  bsub -q normal -n 1  -R "span[hosts=1] rusage[mem=5000]" -M 6000 -o samtools_bt2.log -e samtools_bt2.err "
    cat $cur_contig_list \
    | sed 's/\\n/ /g' \
    | xargs -d '\n' -n 1 -P 40 -I {} /bin/bash -c ' \
        bcftools mpileup -q \$2 -d 800 --threads 1 -f \$3 -r \$1 \${4} > \${1}_mq\${2}.bt2.pileup \
      ' -- {} 1 $refcur ${cwd}/${sample}/${sample}.sorted.RG.bt2.bam

    bcftools concat --threads 38 *_mq1.bt2.pileup \
    | bcftools call --threads 38 -m -v -Ov --ploidy-file ${ploidy} \
    > ${sample}_mq1.bt2.vcf
    vcftools --vcf ${sample}_mq1.bt2.vcf --keep-only-indels --recode --recode-INFO-all --out ${sample}_mq1_onlyindels.bt2.vcf

    vcf2bed --do-not-sort --insertions <${sample}_mq1_onlyindels.bt2.vcf.recode.vcf> ${sample}_mq1_onlyindels.bt2.vcf.recode.insertions.bed
    vcf2bed --do-not-sort --deletions <${sample}_mq1_onlyindels.bt2.vcf.recode.vcf> ${sample}_mq1_onlyindels.bt2.vcf.recode.deletions.bed
    "
done

# GATK HaplotypeCaller
regions='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.ranges.formatted.txt'
jv="--java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=1'"
#jv="--java-options -Xms2G -Xmx2G -XX:ParallelGCThreads=1"
cd $cwd
for sample in ${samples[@]}; do
  cd $cwd/$sample
  cd indels
  mkdir gatk_hc
  cd gatk_hc
  mkdir mm2
  mkdir bt2

#  # Step 1
#  while read r; do
#    r2=$( echo $r | sed 's/:/_/g' | sed 's/-/_/g'  )
#    cd mm2
#    bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo ${r2}.log -eo ${r2}.err "
#      gatk $jv \
#        HaplotypeCaller \
#        -R $refcur \
#        -I ${cwd}/${sample}/${sample}.sorted.RG.bam \
#        -O ${r2}.vcf.gz \
#        --native-pair-hmm-threads 1 \
#        -L $r
#    "
#    cd ..
#
#    cd bt2
#    bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo ${r2}.log -eo ${r2}.err "
#      gatk $jv \
#        HaplotypeCaller \
#        -R $refcur \
#        -I ${cwd}/${sample}/${sample}.sorted.RG.bt2.bam \
#        -O ${r2}.vcf.gz \
#        --native-pair-hmm-threads 1 \
#        -L $r
#    "
#    cd ..
#  done < $regions

  # Step 2
  for dir in mm2 bt2; do
    cd $dir
    if [[ $dir == mm2 ]]; then
      mq_value='50.00'
    else
      mq_value='35.00'
    fi
    bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo indel.log -eo indel.err "
        ls CUR*g.vcf.gz > variants.list
        gatk $jv CombineGVCFs \
          -R $refcur \
          -V variants.list  \
          -O TMP.g.vcf.gz
        gatk $jv GenotypeGVCFs \
          -R $refcur \
          -V TMP.g.vcf.gz \
          -O TMP.vcf
        gatk $jv SelectVariants \
          -R $refcur \
          -V TMP.vcf \
          --select-type-to-include INDEL \
          -O TMP.indel.vcf
        gatk $jv VariantFiltration \
         -R $refcur \
         -V TMP.indel.vcf \
         -O TMP.indel.filter.vcf \
         --filter-expression \"QUAL < 0 || QD < 2.00 || FS > 60.000 || SOR > 3.000  || MQ < $mq_value || MQRankSum < -10.00 || ReadPosRankSum < -6.000 || ReadPosRankSum > 4.000\" \
         --filter-name \"indel_filter\" \
         --filter-expression \"DP < 50 || DP > 300\" \
         --filter-name \"DP_filter\"
        grep -E '^#|PASS' TMP.indel.filter.vcf > indels.vcf
        vcf2bed --do-not-sort --insertions <indels.vcf> insertions.bed
        vcf2bed --do-not-sort --deletions <indels.vcf> deletions.bed
#        rm TMP.*
    "
    cd ..
  done
done


# Merge all the called indel calls. This will give us all regions where identified
# sSNVs would be of low confidence because of alignment errors.
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
samples=( "MUT_11_1" "MUT_15" "WT_1" "WT_19" )
for sample in ${samples[@]}; do
  echo $sample
  cd $cwd/$sample/indels
  rm insertions.list deletions.list

  ls manta/*/*insertions.bed >> insertions.list
  ls samtools/*/*insertions.bed >> insertions.list
  ls gatk_hc/*/*insertions.bed >> insertions.list
  while read r; do
    cat $r >> insertions.all.bed
  done < insertions.list
  cut -f1,2,3 insertions.all.bed > insertions.bed
  sortBed -faidx /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/variantCorrected_AND_kmer_masked_d100-kmer80bp-8_4_Currot_on_Original_v1.0.fasta.fai -i insertions.bed > insertions.sorted.bed
  mergeBed -i insertions.sorted.bed > insertions.sorted.merged.bed

  ls manta/*/*deletions.bed >> deletions.list
  ls samtools/*/*deletions.bed >> deletions.list
  ls gatk_hc/*/*deletions.bed >> deletions.list
  while read r; do
    cat $r >> deletions.all.bed
  done < deletions.list
  cut -f1,2,3 deletions.all.bed > deletions.bed
  sortBed -faidx /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/variantCorrected_AND_kmer_masked_d100-kmer80bp-8_4_Currot_on_Original_v1.0.fasta.fai -i deletions.bed > deletions.sorted.bed
  mergeBed -i deletions.sorted.bed > deletions.sorted.merged.bed

  ## Add buffer of 20 BP around indel positions
  cat insertions.sorted.merged.bed \
  | awk '{printf $1"\t"} {if($2<20) {printf 0"\t"} else {printf $2-20"\t"}} {print $3+20}' \
  > insertions.sorted.merged.padded_20.bed
  cat deletions.sorted.merged.bed \
  | awk '{printf $1"\t"} {if($2<20) {printf 0"\t"} else {printf $2-20"\t"}} {print $3+20}' \
  > deletions.sorted.merged.padded_20.bed

  rm *.all.bed *.sorted.bed insertions.bed deletions.bed
done

## STEP 5a : STRELKA BASED sSNVs
# Calling variants using Strelka: For each sample, call
# variant in that sample (considering it as tumor) compared
# to all other samples (considered as normals). Variants
# identified in the sample in all comparisons are selected
# for further processing. This will be repeated for both
# minimap2 and bowtie2 based alignments.

indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/strelka/'
stlka_config='/srv/netscratch/dep_mercier/grp_schneeberger/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py'
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
samples=( "MUT_11_1" "MUT_15" "WT_1" "WT_19" )

## Run strelka on all combinations of samples
for sample in ${samples[@]}; do
  cd $cwd
  mkdir $sample
  cd $sample
  echo "IN $sample"
  for sample2 in ${samples[@]}; do
    if [[ $sample == $sample2 ]]; then
      continue
    else
      ## Submit MM2 based jobs
      bsub -q multicore40 -n 20 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo ${sample}_${sample2}_mm2.log -eo ${sample}_${sample2}_mm2.err "
        /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/py2.6mg/bin/python \
          $stlka_config \
            --tumourBam ${indir}/${sample}/${sample}.sorted.bam \
            --normalBam ${indir}/${sample2}/${sample2}.sorted.bam \
            --referenceFasta $refcur \
            --runDir ${sample}_vs_${sample2}_mm2 \
            --scanSizeMb 2
        cd ${sample}_vs_${sample2}_mm2
        /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/py2.6mg/bin/python \
          ./runWorkflow.py \
            -m local \
            -j 40 \
            -g 50
     "

      ## Submit BT2 based jobs
      bsub -q multicore40 -n 20 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo ${sample}_${sample2}_bt2.log -eo ${sample}_${sample2}_bt2.err "
        /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/py2.6mg/bin/python \
          $stlka_config \
            --tumourBam ${indir}/${sample}/${sample}.sorted.bt2.bam \
            --normalBam ${indir}/${sample2}/${sample2}.sorted.bt2.bam \
            --referenceFasta $refcur \
            --runDir ${sample}_vs_${sample2}_bt2 \
            --scanSizeMb 2
        cd ${sample}_vs_${sample2}_bt2
        /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/py2.6mg/bin/python \
          ./runWorkflow.py \
            -m local \
            -j 40 \
            -g 50
      "
    fi
  done
done



## STEP 5b : Mutect2 BASED sSNVs
# First call mutations in all individual samples, RojoPassion,
# and the parents data in tumor-only mode. Then for each sample
# create a Panel_of_Normals using all vcf files (other than the
# focal sample) identified above. Finally, call variants in the
# comparative mode using sample data (from cell reads) and the
# corresponding Panel_of_Normals.
#---------------- ADD RG TAGS TO BAM FILES ---------------------------------
# Add RG tag in the bam files for the four scDNAseq samples
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/mutect2/'
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
picard='/srv/netscratch/dep_mercier/grp_schneeberger/software/picard-tools-1.131/picard.jar'
regions='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.ranges.formatted.txt'

samples=( "MUT_11_1" "MUT_15" "WT_1" "WT_19" )
cd $cwd
for sample in ${samples[@]}; do
  cd $cwd
  mkdir $sample
  cd $sample
    bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 8000 -oo ${sample}_add_RG_mm2.log -eo ${sample}_${add_RG}_mm2.err "
      java -XX:ParallelGCThreads=1 -jar $picard AddOrReplaceReadGroups I=${indir}/${sample}/${sample}.sorted.bam O=${indir}/${sample}/${sample}.sorted.RG.bam LB=${sample} PL=illumina PU=unit1 SM=${sample}
      samtools index ${indir}/${sample}/${sample}.sorted.RG.bam
    "

    bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 8000 -oo ${sample}_add_RG_bt2.log -eo ${sample}_${add_RG}_bt2.err "
      java -XX:ParallelGCThreads=1 -jar $picard AddOrReplaceReadGroups I=${indir}/${sample}/${sample}.sorted.bt2.bam O=${indir}/${sample}/${sample}.sorted.RG.bt2.bam LB=${sample} PL=illumina PU=unit1 SM=${sample}
      samtools index ${indir}/${sample}/${sample}.sorted.RG.bt2.bam
    "
done


# Add RG tag in the bam files for the Currot and OrangeRed samples
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/'
samples=( 'currot' 'orangered' )
for sample in ${samples[@]}; do
  cd ${cwd}/$sample
  bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 8000 -oo ${sample}_add_RG_mm2.log -eo ${sample}_${add_RG}_mm2.err "
      java -XX:ParallelGCThreads=1 -jar $picard AddOrReplaceReadGroups \
        I=${sample}_cur_v1.1.dedup.sorted.bam \
        O=${sample}_cur_v1.1.dedup.sorted.RG.bam \
        LB=${sample} PL=illumina PU=unit1 SM=${sample}
      samtools index ${sample}_cur_v1.1.dedup.sorted.RG.bam
    "
done

# Add RG tag in the bam files for Rojo Passion merged samples also submit jobs for calling variants
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
regions='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.ranges.formatted.txt'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/rp_leaf_illumina/merged_samples/'
picard='/srv/netscratch/dep_mercier/grp_schneeberger/software/picard-tools-1.131/picard.jar'
cd $cwd
# Step 1
bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 8000 -oo add_RG.log -eo add_RG.err "
    java -XX:ParallelGCThreads=1 -jar $picard AddOrReplaceReadGroups \
      I=merged.bam \
      O=merged.RG.bam \
      LB=RP PL=illumina PU=unit1 SM=RP
    samtools index merged.RG.bam
    while read r; do
      r2=\$( echo \$r | sed 's/:/_/g' | sed 's/-/_/g'  )
      bsub -q multicore20 -n 1 -R \"span[hosts=1] rusage[mem=5000]\" -M 5000 -oo \${r2}_mm2.log -eo \${r2}_mm2.err \"
        gatk --java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=1' \
          Mutect2 \
            -R $refcur \
            -I merged.RG.bam \
            -tumor RP \
            -O \${r2}.vcf.gz \
            --native-pair-hmm-threads 1 \
            -L \$r
      \"
  done < $regions
"
# Step 2
bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo sort.log -eo sort.err "
  bcftools concat -O z CUR*.vcf.gz > all_variants.vcf.gz
  java -jar $picard SortVcf I=all_variants.vcf.gz O=all_variants.sorted.vcf.gz
"

#---------------- DETECT VARIANTS ---------------------------------
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/mutect2/'
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
regions='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.ranges.formatted.txt'
# Get variations in tumor-only mode for scdna-seq samples
cd $cwd
for sample in ${samples[@]}; do
  cd $cwd
  mkdir $sample
  cd $sample
  mkdir tumor_only_mode
  cd tumor_only_mode

  # Step 1
  while read r; do
    r2=$( echo $r | sed 's/:/_/g' | sed 's/-/_/g'  )
    bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo ${r2}_mm2.log -eo ${r2}_mm2.err "
      gatk --java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=1' \
        Mutect2 \
        -R $refcur \
        -I ${indir}/${sample}/${sample}.sorted.RG.bam \
        -tumor $sample \
        -O ${r2}_mm2.vcf.gz \
        --native-pair-hmm-threads 1 \
        -L $r
    "

     bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo ${r2}_bt2.log -eo ${r2}_bt2.err "
      gatk --java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=1' \
        Mutect2 \
        -R $refcur \
        -I ${indir}/${sample}/${sample}.sorted.RG.bt2.bam \
        -tumor $sample \
        -O ${r2}_bt2.vcf.gz \
        --native-pair-hmm-threads 1 \
        -L $r
    "
  done < $regions

  # Step 2
  bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo sort.log -eo sort.err "
    bcftools concat -O z CUR*_mm2.vcf.gz > all_variants.mm2.vcf.gz
    java -jar $picard SortVcf I=all_variants.mm2.vcf.gz O=all_variants.mm2.sorted.vcf.gz

    bcftools concat -O z CUR*_bt2.vcf.gz > all_variants.bt2.vcf.gz
    java -jar $picard SortVcf I=all_variants.bt2.vcf.gz O=all_variants.bt2.sorted.vcf.gz
  "
done


# Get variations in tumor-only mode for the Currot and OrangeRed samples
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/'
samples=( 'currot' 'orangered' )
for sample in ${samples[@]}; do
  cd ${cwd}/$sample
  mkdir mutect2_var_call
  cd mutect2_var_call

  # Step 1
  while read r; do
    r2=$( echo $r | sed 's/:/_/g' | sed 's/-/_/g'  )
    bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo ${r2}_mm2.log -eo ${r2}_mm2.err "
      gatk --java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=1' \
        Mutect2 \
          -R $refcur \
          -I ${cwd}/${sample}/${sample}_cur_v1.1.dedup.sorted.RG.bam \
          -tumor $sample \
          -O ${sample}.${r2}.vcf.gz \
          --native-pair-hmm-threads 1 \
          -L $r
    "
  done < $regions

  # Step 2
  bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo sort.log -eo sort.err "
    bcftools concat -O z ${sample}.CUR*.vcf.gz > all_variants.vcf.gz
    java -jar $picard SortVcf I=all_variants.vcf.gz O=all_variants.sorted.vcf.gz
  "
done

# Call variations using Mutect2
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/mutect2/'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
regions='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.ranges.formatted.txt'
samples=( "MUT_11_1" "MUT_15" "WT_1" "WT_19" )
cd $cwd

rf  */tumor_only_mode/all_variants.mm2.sorted.vcf.gz > variants_mm2.args
rf  */tumor_only_mode/all_variants.bt2.sorted.vcf.gz > variants_bt2.args
echo /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/rp_leaf_illumina/merged_samples/all_variants.sorted.vcf.gz | tee -a variants_mm2.args variants_bt2.args
echo /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/mutect2_var_call/all_variants.sorted.vcf.gz | tee -a variants_mm2.args variants_bt2.args
echo /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/orangered/mutect2_var_call/all_variants.sorted.vcf.gz | tee -a variants_mm2.args variants_bt2.args

# Create Panel Of Normals
for sample in ${samples[@]}; do
  cd $cwd
  cd $sample
  for dir in mm2 bt2; do
    mkdir $dir
    cd $dir
    bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo pon.log -eo pon.err "
      gatk --java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=1 -XX:ActiveProcessorCount=1' \
        CreateSomaticPanelOfNormals \
        -vcfs ${cwd}/variants_${dir}.args \
        -O pon_${dir}.vcf.gz
    "
    cd ..
  done
done

# Get variants
for sample in ${samples[@]}; do
  cd $cwd
  cd $sample
  for dir in mm2 bt2; do
    echo $sample $dir
    tool=''
    if [[ $dir == 'bt2' ]]; then
      tool='.bt2'
    fi
    cd $dir
    for sample2 in ${samples[@]}; do
      if [[ $sample == $sample2 ]]; then
        continue
      else

        mkdir ${sample}_vs_${sample2}
        cd ${sample}_vs_${sample2}
        while read r; do
          if [[ $r != "CUR"* ]]; then
            continue
          else
            r2=$( echo $r | sed 's/:/_/g' | sed 's/-/_/g'  )
            bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=5000]" -M 5000 -oo ${r2}.log -eo ${r2}.err "
              gatk --java-options '-Xms5G -Xmx5G -XX:ParallelGCThreads=1 -XX:ActiveProcessorCount=1' \
                Mutect2 \
                  -R $refcur \
                  -I ${indir}/${sample}/${sample}.sorted.RG${tool}.bam \
                  -tumor $sample \
                  -I ${indir}/${sample2}/${sample2}.sorted.RG${tool}.bam \
                  -normal ${sample2} \
                  -O ${r2}.vcf.gz \
                  --panel-of-normals ${cwd}/${sample}/${dir}/pon_${dir}.vcf.gz \
                  --native-pair-hmm-threads 1 \
                  -L $r
            "
          fi
        done < $regions
        cd ..
      fi
    done
    cd ..
  done
done


#######################################################################################
######################## END OF NECESSARY COMMANDS ####################################
#######################################################################################
## Testing commands
cd $cwd
samtools depth -a -q 30 -Q 10 -d 0 WT_1/WT_1.sorted.bam WT_19/WT_19.sorted.bam MUT_11_1/MUT_11_1.sorted.bam MUT_15/MUT_15.sorted.bam > read_depth_q30_Q10.txt

#-----------------------------------------------------------------------
## Test different sSNP identification methods using simulated BAM files
picard='/srv/netscratch/dep_mercier/grp_schneeberger/software/picard-tools-1.131/picard.jar'
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/CUR1G.fasta'
curidxbt2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/CUR1G'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/'
cd $cwd

bowtie2 --end-to-end \
        --very-sensitive \
        --threads 40 \
        -x $curidxbt2 \
        -1 WT_1_R1.fastq.gz \
        -2 WT_1_R2.fastq.gz \
        -S CUR1G.bt2.sam

samtools sort -O BAM CUR1G.bt2.sam > CUR1G.bt2.sorted.bam
samtools index CUR1G.bt2.sorted.bam
samtools view -b -f 2 CUR1G.bt2.sorted.bam > CUR1G.bt2.sorted.onlymapped.bam

# Add RG information
java -jar $picard AddOrReplaceReadGroups I=CUR1G.bt2.sorted.onlymapped.bam O=CUR1G.bt2.sorted.onlymapped.RG.bam LB=WT_1 PL=illumina PU=unit1 SM=1

# Check validity of BAM file
java -jar $picard ValidateSamFile I=CUR1G.bt2.sorted.onlymapped.RG.bam O=validation.report M=VERBOSE

samtools index CUR1G.bt2.sorted.onlymapped.RG.bam

# Get all positions with no alternate base and read-coverage between 100-200 bps, SNPs would be added within these positions only
bam-readcount -w 0 -f $refcur CUR1G.bt2.sorted.onlymapped.RG.bam \
| awk '{if($4>=100 && $4<=200) {n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":");  print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2], f[1], f[2], g[1], g[2], h[1], h[2]}}' \
 > bam_read_counts_b30_q20.bt2.txt

python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/get_conserved_positions.py bam_read_counts_b30_q20.bt2.txt

addsnv.py -v mutate_loci_bam_read_counts_b30_q20.bt2.txt \
  -f CUR1G.bt2.sorted.onlymapped.RG.bam \
  -r ../CUR1G.fasta \
  -o CUR1G.bt2.sorted.onlymapped.RG.spiked.bam \
  -p 60 \
  --picardjar /srv/netscratch/dep_mercier/grp_schneeberger/software/picard-tools-1.131/picard.jar \
  --mindepth 100 \
  --maxdepth 200 \
  --tagreads \
  --aligner bowtie2 \
  --alignopts end-to-end:,very-sensitive:,threads:40,bowtie2ref:CUR1G \
  --seed 1

cd WT_1
samtools sort -O BAM -@ 50 CUR1G.bt2.sorted.onlymapped.RG.spiked.bam > CUR1G.spiked.sorted.bam
samtools index CUR1G.spiked.sorted.bam

## ABOVE STEPS ARE REPEATED FOR WT_19 SAMPLE AS WELL

# Using py2.6mg environment run STRELKA
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/'
cd $cwd
python /srv/netscratch/dep_mercier/grp_schneeberger/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py \
 --normalBam WT_19/CUR1G.bt2.sorted.onlymapped.RG.bam \
 --tumourBam WT_1/CUR1G.spiked.sorted.bam \
 --referenceFasta CUR1G.fasta \
 --runDir ./mut_in_wt_1 \
 --scanSizeMb 2
cd mut_in_wt_1
python ./runWorkflow.py \
  -m local \
  -j 60 \
  -g 250

# run strelka for simulated WT_19
python /srv/netscratch/dep_mercier/grp_schneeberger/software/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --tumourBam CUR1G.spiked.sorted.bam --normalBam ../WT_1/CUR1G.bt2.sorted.onlymapped.RG.bam --referenceFasta ../CUR1G.fasta --runDir ./mut_in_wt_19 --scanSizeMb 2
python ./runWorkflow.py \
  -m local \
  -j 60 \
  -g 250

# Running mutect2 on WT_1
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/WT_1
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
gatk Mutect2 \
  -R $refcur \
  -I CUR1G.spiked.sorted.RG.bam \
  -I ../WT_19/CUR1G.spiked.sorted.RG.bam \
  -tumor-sample WT_1 \
  -normal WT_19 \
  --native-pair-hmm-threads 40 \
  -O WT_1_vs_WT_19.vcf.gz

# Running mutect2 on WT_19
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/tool_based_analysis/test_on_simulated_data/WT_19
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
gatk Mutect2 \
  -R $refcur \
  -I CUR1G.spiked.sorted.RG.bam \
  -I ../WT_1/CUR1G.spiked.sorted.RG.bam \
  -tumor-sample WT_19 \
  -normal WT_1 \
  --native-pair-hmm-threads 40 \
  -O WT_19_vs_WT_1.vcf.gz
#-----------------------------------------------------------------------
#######################################################################################
### TEMPORARY COMMANDS (FOR RUNNING ON DELL-NODES OR RUNNING SPECIFIC SUB-COMMANDS ####
#######################################################################################
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_1/
# Get pileup data


gatk VariantsToTable \
  -R $refcur \
  -V TMP.indel.vcf \
  -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F MQRankSum -F FS -F ReadPosRankSum -F SOR \
  -O TMP.indel.table


gatk CombineGVCFs \
          -R $refcur \
          -V variants.list  \
          -O TMP.g.vcf.gz
        gatk GenotypeGVCFs \
          -R $refcur \
          -V TMP.g.vcf.gz \
          -O TMP.vcf
        gatk SelectVariants \
          -R $refcur \
          -V TMP.vcf \
          --select-type-to-include INDEL \
          -O TMP.indel.vcf
        gatk VariantFiltration \
         -R $refcur \
         -V TMP.indel.vcf \
         -O TMP.indel.filter.vcf \
         --filter-expression \"QUAL < 0 || QD < 2.00 || FS > 60.000 || SOR > 3.000  || MQ < $mq_value || MQRankSum < -10.00 || ReadPosRankSum < -6.000 || ReadPosRankSum > 4.000\" \
         --filter-name \"indel_filter\" \
         --filter-expression \"DP < 50 || DP > 300\" \
         --filter-name \"DP_filter\"
        grep -E '^#|PASS' TMP.indel.filter.vcf > indels.vcf
        vcf2bed --do-not-sort --insertions <indels.vcf> insertions.bed
        vcf2bed --do-not-sort --deletions <indels.vcf> deletions.bed


gatk VariantFiltration \
  -R $refcur \
  -V TMP.indel.vcf \
  -O TMP.indel.filter.vcf \
  --filter-expression "QUAL < 0 || QD < 2.00 || FS > 60.000 || SOR > 3.000  || MQ < $mq_value || MQRankSum < -10.0 || ReadPosRankSum < -6.000 || ReadPosRankSum > 4.000" \
  --filter-name "indel_filter" \
  --filter-expression "DP < 50 || DP > 300" \
  --filter-name "DP_filter"
        grep -E '^#|PASS' TMP.indel.filter.vcf > indels.vcf
        vcf2bed --do-not-sort --insertions <indels.vcf> insertions.bed
        vcf2bed --do-not-sort --deletions <indels.vcf> deletions.bed