## Step 1: Index current Currot and OrangeRed assemblies to be used as reference genomes for barcode correction

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
    bsub -q multicore40 -n 20 -R "span[hosts=1] rusage[mem=80000]" -M 80000 -oo align_cells.log -eo align_cells.err "
    xargs -a barcodes_list -n 1 -P 20 -I {} /bin/bash -c '
        cd {}
        bc=\$(basename {})
#        minimap2 -ax sr -t 1 \$2 \${bc}_R1.fastq.gz \${bc}_R2.fastq.gz \
#        | samtools view -O BAM - \
#        | samtools sort - > \${bc}.sorted.bam
#        samtools index \${bc}.sorted.bam

        ## Remove duplicated reads
        samtools sort -n -O BAM \${bc}.sorted.bam \
        | samtools fixmate -c -m -O BAM - - \
        | samtools sort -O BAM - \
        | samtools markdup -S -s -O BAM - - \
        | samtools sort -O BAM - > \${bc}.DUPmarked.bam
        samtools index \${bc}.DUPmarked.bam

        # Get non-dup reads in fastq.gz from the markdup bam files
        python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/remove_duplicates_from_marked_bam_file.py \${bc}.DUPmarked.bam -p \${bc}_dedup

        # Get non-dup bam file from the markdup bam file
        samtools view -G 1024 -O BAM \${bc}.DUPmarked.bam > \${bc}.DUPmarked.deduped.bam
        samtools index \${bc}.DUPmarked.deduped.bam

        # Get BAM coverage
        bamCoverage --numberOfProcessors 1 -b \${bc}.DUPmarked.deduped.bam -of bedgraph -o \${bc}.bedgraph ' -- {} $2
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

#############################################################################################################
############################ Step  3: Variant Calling from pooled data  #####################################
#############################################################################################################
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/get_cells/all_barcodes/'
curidx='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta.mm2_Xsr.idx'
samples=( "MUT_11_1" "MUT_15" "WT_1" "WT_19" )

cd $cwd

## Step 3a: Merge individual barcodes fastq to one sample fastq
for sample in ${samples[@]}; do
  cd $cwd
  mkdir $sample
  cd $sample
  while read r; do
    bc=$(basename $r)
    cat ${indir}${sample}/barcodes/${bc}/${bc}_dedup_R1.fastq.gz >> ${sample}_R1.fastq.gz
    cat ${indir}${sample}/barcodes/${bc}/${bc}_dedup_R2.fastq.gz >> ${sample}_R2.fastq.gz
  done < ${indir}${sample}/barcodes_list
done

## Step 3b: Align merged reads against the currot reference genome
for sample in ${samples[@]}; do
  cd $cwd
  mkdir $sample
  cd $sample
  bsub -q multicore40 -n 20 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo align_reads_to_cur.log -eo align_reads_to_cur.err "
    minimap2 -ax sr -t 18 $curidx ${sample}_R1.fastq.gz ${sample}_R2.fastq.gz \
    | samtools view -O BAM - \
    | samtools sort - > ${sample}.sorted.bam
    samtools index ${sample}.sorted.bam

    # Duplicates were removed from the indivual cells. Therefore, there is no need to remove duplicates from the merged data
#    samtools sort -n -@ 18 -O BAM ${sample}.sorted.bam \
#    | samtools fixmate -c -m -O BAM -@ 18 - - \
#    | samtools sort -O BAM -@ 18 - \
#    | samtools markdup -r -S -s -O BAM -@ 18 - - \
#    | samtools sort -O BAM -@ 18 - > ${sample}.rmdup.sorted.bam
#    samtools index ${sample}.rmdup.sorted.bam
  "
done

## Testing commands
cd $cwd
samtools depth -a -q 30 -Q 10 -d 0 WT_1/WT_1.rmdup.sorted.bam WT_19/WT_19.rmdup.sorted.bam MUT_11_1/MUT_11_1.rmdup.sorted.bam MUT_15/MUT_15.rmdup.sorted.bam > read_depth_q30_Q10.txt

## Step 3c: getting variant positions
refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta'
test='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/WT_1/head.bam'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/'
samples=( "MUT_11_1" "MUT_15" "WT_1" "WT_19" )


for sample in ${samples[@]}; do
  cd $cwd
  cd $sample
# GET READCOUNT FROM HIGH QUALITY READS AND BASES
#  bam-readcount -b 30 -q 20 -w 0 -f $refcur *.rmdup.sorted.bam| awk '{if($4>3) {n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":"); n6=split($11,f, ":"); n7=split($12,g, ":"); n8=split($13,h, ":");  print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2], f[1], f[2], g[1], g[2], h[1], h[2]}}' > bam_read_counts_b30_q20.txt &

  # GET POSITIONS WITH AT LEAST THREE NON-REFERENCE BASES
  python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py bam_read_counts_b30_q20.txt &
done


#bam-readcount -b 30 -q 20 -w 0 -f $refcur $test |  awk '{if($4>3) {n1=split($6,a,":"); n2=split($7,b,":"); n3=split($8,c, ":"); n4=split($9,d,":"); n5=split($10,e,":"); n6=split($11,f, ":"); n7=split($12,g, ":"); n8=split($13,h, ":");  print $1, $2, $3, $4, a[2], b[2], c[2], d[2], e[2], f[1], f[2], g[1], g[2], h[1], h[2]}}'
#
#awk '{if($4>3) n1=split(\$6,a,\":\"); n2=split(\$7,b,\":\"); n3=split(\$8,c, \":\"); n4=split(\$9,d,\":\"); n5=split(\$10,e,\":\"); n6=split(\$11,f, \":\"); print \$1, \$2, \$3, \$4, a[1], a[2], b[1], b[2], c[1], c[2], d[1], d[2], e[1], e[2], f[1], f[2]}'




samtools sort -n -O BAM ${bc}.sorted.bam \
        | samtools fixmate -c -m -O BAM - - \
        | samtools sort -O BAM - \
        | samtools markdup -r -S -s -O BAM - - \
        | samtools sort -O BAM - > ${bc}.rmdup.sorted.bam


