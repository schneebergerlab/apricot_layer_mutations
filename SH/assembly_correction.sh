########################### Polish the assemblies ##############################

## CURROT assembly

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/cur/
curR1='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/currot_ql-trimmed-pair1.fastq.gz'
curR2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/currot_ql-trimmed-pair2.fastq.gz'
picard=/srv/netscratch/dep_mercier/grp_schneeberger/bin/gatk/picard.jar

# Build-indices
bowtie2-build cur.v3.fasta cur.v3.fasta &> bowtie2_index.log&
bwa index -a bwtsw cur.v3.fasta &>bwa_index.log&
samtools faidx cur.v3.fasta &>samtools_index.log&

# Align parent reads, mark duplicates and get insert_size
bsub -q ioheavy -n 40 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo ill_aln.log -eo ill_aln.err "
  bowtie2 -x cur.v3.fasta -1 ${curR1} -2 ${curR2} -p 40 \
  | samtools sort -@ 40 -O BAM - \
  > cur.v3.cursr.sorted.bam
  samtools index -@ 40 cur.v3.cursr.sorted.bam
  java -Xmx45G -XX:ParallelGCThreads=1 -XX:ActiveProcessorCount=40 \
    -jar $picard MarkDuplicates \
    I=cur.v3.cursr.sorted.bam \
    O=cur.v3.cursr.sorted.markdup.bam \
    M=cur.v3.cursr.sorted.markdup.stats.txt
  java -Xmx45G -XX:ParallelGCThreads=1 -XX:ActiveProcessorCount=40 \
    -jar $picard CollectInsertSizeMetrics \
    I=cur.v3.cursr.sorted.markdup.bam \
    O=cur.v3.cursr.sorted.markdup.insert_size.txt \
    H=cur.v3.cursr.sorted.markdup.insert_size.pdf
  samtools index cur.v3.cursr.sorted.markdup.bam
"

# Get contig alignments and fasta sequence
cat cur.v3.fasta.chrsize \
| cut -f1 \
| xargs -n 1 -P 60 -I {} bash -c '
  mkdir -p contigs/${1}
#  samtools view -b cur.v3.cursr.sorted.markdup.bam ${1} > contigs/${1}/${1}.sorted.markdup.bam
#  samtools index contigs/${1}/${1}.sorted.markdup.bam
  hometools getchr --chrs ${1} -o contigs/${1}/${1}.fasta cur.v3.fasta
' -- {}


# Run pilon for polishing using illumina reads
pilon=/opt/share/software/packages/pilon-1.22/pilon.jar
contigs=$(cat cur.v3.fasta.chrsize | cut -f1)
for contig in ${contigs[@]}; do
  cd contigs/${contig}
  bsub -q short -n 3 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo pilon.log -eo pilon.err "
    java -Xmx9G -XX:ActiveProcessorCount=3 \
      -jar $pilon \
              --genome ${contig}.fasta \
              --fix bases \
              --frags ${contig}.sorted.markdup.bam \
              --threads 4 \
              --outdir ./ \
              --output pilon_${contig} \
              --changes \
              --mindepth 0.8 \
              --minmq 1 \
              --minqual 10 \
              --chunksize 15000000
  "
  cd ../../
done

# Merge polished contigs
cat contigs/*/pilon_*.fasta >> cur.v3.pilon_polished.fasta
hometools seqsize cur.v3.pilon_polished.fasta > cur.v3.pilon_polished.fasta.chrsize

# Align HiFi-reads
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/'
minimap2 -ax map-pb -t 60 --secondary=no \
  cur.v3.pilon_polished.fasta \
  ${indir}/haplotype-cur.fasta.gz \
  ${indir}/haplotype-unknown.fasta.gz \
| samtools sort -@ 60 -O BAM - \
> cur.v3.pilon_polished.sorted.bam
samtools index cur.v3.pilon_polished.sorted.bam

# Extract contig alignments and reads. NOTE: racon do not accept BAM and requires SAM as input
cat cur.v3.fasta.chrsize \
| cut -f1 \
| xargs -n 1 -P 60 -I {} bash -c '
    echo ${1}
    samtools view -h cur.v3.pilon_polished.sorted.bam ${1}_pilon > contigs/${1}/${1}.hifi.sorted.sam
    samtools fasta contigs/${1}/${1}.hifi.sorted.sam | gzip > contigs/${1}/${1}.hifi.fasta.gz
    gzip contigs/${1}/${1}.hifi.sorted.sam
' -- {}



# Polish using Racon and hifi-reads:
# racon outputs fasta id with id, e.g., "homLG_6_LG_26_utg000001l_pilon LN:i:1135064 RC:i:2291 XC:f:0.999560"
# LN - sequence length after polishing
# RC - number of reads used for polishing the sequence
# XC - percentage of polished windows in the sequnce
# why no trimming option: --no-trimming, see here: https://github.com/isovic/racon/issues/126
# 20201026: better use trim as it is that at ends of contigs, there is always low coverage according to short read alignment.
contigs=$(cat cur.v3.fasta.chrsize | cut -f1)
for contig in ${contigs[@]}; do
  cd contigs/${contig}
  if [ $contig]
  bsub -q short -n 3 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo rocon.log -eo racon.err "
    racon -u --no-trimming -t 3 \
      ${contig}.hifi.fasta.gz \
      ${contig}.hifi.sorted.sam.gz \
      pilon_${contig}.fasta > racon_${contig}.no_trimming.fasta.TMP
  "
  cd ../../
done

cat contigs/*/racon_*.no_trimming.fasta > cur.v3.racon_polished.fasta
hometools seqsize cur.v3.racon_polished.fasta > cur.v3.racon_polished.fasta.chrsize

# check size change: lg=; fasta_length ../polished_ref_48_LGs_id_updated/polish_${lg}/illu_polished_${lg}.fasta | grep  '>' | awk '{s+=$2} END {print s}'; fasta_length ./racon_polish_${lg}/illu_hifi_polished_${lg}.fasta | grep  '>' | awk '{s+=$5} END {print s}'

>zracon_polished_size_versus_original_size
while read lg; do
  echo "illu_"${lg} >> zracon_polished_size_versus_original_size
  fasta_length ../polished_ref_48_LGs_id_updated/polish_${lg}/illu_polished_${lg}.fasta | grep  '>' | awk '{s+=$2} END {print s}' >> zracon_polished_size_versus_original_size
  echo "illu_plus_hifi_"${lg} >> zracon_polished_size_versus_original_size
  fasta_length ./racon_polish_${lg}/illu_hifi_polished_${lg}.fasta | grep  '>' | awk '{s+=$5} END {print s}' >> zracon_polished_size_versus_original_size
done < ../../s7_selected_long_contigs_LG_wise_assembly/asm_version_20200910/fa_to_run.list

# total size after illu_hifi_polished: 3066491438





## ORANGERED assembly

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/ora/
oraR1='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/orangered/orangered_ql-trimmed-pair1.fastq.gz'
oraR2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/orangered/orangered_ql-trimmed-pair2.fastq.gz'
picard=/srv/netscratch/dep_mercier/grp_schneeberger/bin/gatk/picard.jar

# Build-indices
bowtie2-build ora.v3.fasta ora.v3.fasta &> bowtie2_index.log&
bwa index -a bwtsw ora.v3.fasta &>bwa_index.log&
samtools faidx ora.v3.fasta &>samtools_index.log&

# Align parent reads, mark duplicates and get insert_size
bsub -q ioheavy -n 40 -R "span[hosts=1] rusage[mem=50000]" -M 60000 -oo ill_aln.log -eo ill_aln.err "
  bowtie2 -x ora.v3.fasta -1 ${oraR1} -2 ${oraR2} -p 40 \
  | samtools sort -@ 40 -O BAM - \
  > ora.v3.orasr.sorted.bam
  samtools index -@ 20 ora.v3.orasr.sorted.bam
  java -Xmx45G -XX:ParallelGCThreads=1 -XX:ActiveProcessorCount=40 \
    -jar $picard MarkDuplicates \
    I=ora.v3.orasr.sorted.bam \
    O=ora.v3.orasr.sorted.markdup.bam \
    M=ora.v3.orasr.sorted.markdup.stats.txt
  java -Xmx45G -XX:ParallelGCThreads=1 -XX:ActiveProcessorCount=40 \
    -jar $picard CollectInsertSizeMetrics \
    I=ora.v3.orasr.sorted.markdup.bam \
    O=ora.v3.orasr.sorted.markdup.insert_size.txt \
    H=ora.v3.orasr.sorted.markdup.insert_size.pdf
  samtools index ora.v3.orasr.sorted.markdup.bam
"

# Get contig alignments
cat ora.v3.fasta.chrsize \
| cut -f1 \
| xargs -n 1 -P 60 -I {} bash -c '
  mkdir -p contigs/${1}
  samtools view -b ora.v3.orasr.sorted.markdup.bam ${1} > contigs/${1}/${1}.sorted.markdup.bam
  samtools index contigs/${1}/${1}.sorted.markdup.bam
  hometools getchr --chrs ${1} -o contigs/${1}/${1}.fasta ora.v3.fasta
' -- {}


# Run pilon for polishing using illumina reads
pilon=/opt/share/software/packages/pilon-1.22/pilon.jar
contigs=$(cat ora.v3.fasta.chrsize | cut -f1)
for contig in ${contigs[@]}; do
  cd contigs/${contig}
  bsub -q short -n 3 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo pilon.log -eo pilon.err "
    java -Xmx9G -XX:ActiveProcessorCount=3 \
      -jar $pilon \
              --genome ${contig}.fasta \
              --fix bases \
              --frags ${contig}.sorted.markdup.bam \
              --threads 4 \
              --outdir ./ \
              --output pilon_${contig} \
              --changes \
              --mindepth 0.8 \
              --minmq 1 \
              --minqual 10 \
              --chunksize 15000000
  "
  cd ../../
done

# Merge polished contigs
cat contigs/*/pilon_*.fasta >> ora.v3.pilon_polished.fasta
hometools seqsize ora.v3.pilon_polished.fasta > ora.v3.pilon_polished.fasta.chrsize

# Align HiFi-reads
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/'
minimap2 -ax map-pb -t 60 --secondary=no \
  ora.v3.pilon_polished.fasta \
  ${indir}/haplotype-ora.fasta.gz \
  ${indir}/haplotype-unknown.fasta.gz \
| samtools sort -@ 60 -O BAM - \
> ora.v3.pilon_polished.sorted.bam
samtools index ora.v3.pilon_polished.sorted.bam


# Extract contig alignments and reads. NOTE: racon do not accept BAM and requires SAM as input
cat ora.v3.fasta.chrsize \
| cut -f1 \
| xargs -n 1 -P 60 -I {} bash -c '
    echo ${1}
    samtools view -h ora.v3.pilon_polished.sorted.bam ${1}_pilon > contigs/${1}/${1}.hifi.sorted.sam
    samtools fasta contigs/${1}/${1}.hifi.sorted.sam | gzip > contigs/${1}/${1}.hifi.fasta.gz
    gzip contigs/${1}/${1}.hifi.sorted.sam
' -- {}


# Polish using Racon and hifi-reads:
# racon outputs fasta id with id, e.g., "homLG_6_LG_26_utg000001l_pilon LN:i:1135064 RC:i:2291 XC:f:0.999560"
# LN - sequence length after polishing
# RC - number of reads used for polishing the sequence
# XC - percentage of polished windows in the sequnce
# why no trimming option: --no-trimming, see here: https://github.com/isovic/racon/issues/126
# 20201026: better use trim as it is that at ends of contigs, there is always low coverage according to short read alignment.
contigs=$(cat ora.v3.fasta.chrsize | cut -f1)
for contig in ${contigs[@]}; do
  cd contigs/${contig}
  bsub -q short -n 3 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo racon.log -eo racon.err "
    racon -u --no-trimming -t 3 \
      ${contig}.hifi.fasta.gz \
      ${contig}.hifi.sorted.sam.gz \
      pilon_${contig}.fasta > racon_${contig}.no_trimming.fasta
  "
  cd ../../
done

cat contigs/*/racon_*.no_trimming.fasta > ora.v3.racon_polished.fasta
hometools seqsize ora.v3.racon_polished.fasta > ora.v3.racon_polished.fasta.chrsize
