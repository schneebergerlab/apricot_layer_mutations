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
  bsub -q short -n 3 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo rocon.log -eo racon.err "
    racon -u --no-trimming -t 3 \
      ${contig}.hifi.fasta.gz \
      ${contig}.hifi.sorted.sam.gz \
      pilon_${contig}.fasta > racon_${contig}.no_trimming.fasta
  "
  cd ../../
done

cat contigs/*/racon_*.no_trimming.fasta > cur.v3.racon_polished.fasta
hometools seqsize cur.v3.racon_polished.fasta > cur.v3.racon_polished.fasta.chrsize


# Racon second run after removing alignments with mapping quality <10 and removing reads with clippings
# Align HiFi-reads for second RACON run
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/'
minimap2 -ax map-pb -t 60 --secondary=no \
  cur.v3.racon_polished.fasta \
  ${indir}/haplotype-cur.fasta.gz \
  ${indir}/haplotype-unknown.fasta.gz \
| samtools sort -@ 60 -O BAM - \
> cur.v3.racon_polished.sorted.bam
samtools view -h -q10 cur.v3.racon_polished.sorted.bam \
| awk '$6 !~ /H|S/{print}' \
| samtools view -bS - \
> cur.v3.racon_polished.q10.noclips.sorted.bam
samtools index -@60 cur.v3.racon_polished.q10.noclips.sorted.bam

# Extract contig alignments and reads. NOTE: racon do not accept BAM and requires SAM as input
cat cur.v3.fasta.chrsize \
| cut -f1 \
| xargs -n 1 -P 60 -I {} bash -c '
    echo ${1}
    mkdir contigs2/${1}
    samtools view -h cur.v3.racon_polished.q10.noclips.sorted.bam ${1}_pilon > contigs2/${1}/${1}.hifi.sorted.sam
    samtools fasta contigs2/${1}/${1}.hifi.sorted.sam | gzip > contigs2/${1}/${1}.hifi.fasta.gz
    gzip contigs2/${1}/${1}.hifi.sorted.sam
    hometools getchr --chrs ${1}_pilon -o contigs2/${1}/${1}.fasta cur.v3.racon_polished.fasta
' -- {}


contigs=$(cat cur.v3.fasta.chrsize | cut -f1)
for contig in ${contigs[@]}; do
  cd contigs2/${contig}
  bsub -q short -n 3 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo racon.log -eo racon.err "
    racon -u --no-trimming -t 3 \
      ${contig}.hifi.fasta.gz \
      ${contig}.hifi.sorted.sam.gz \
      ${contig}.fasta > racon_${contig}.no_trimming.fasta
  "
  cd ../../
done

contigs=$(cat cur.v3.fasta.chrsize | cut -f1)
rm cur.v3.racon_polished2.fasta
for contig in ${contigs[@]}; do
  lc=$(wc -l contigs2/${contig}/racon_${contig}.no_trimming.fasta | cut -d' ' -f1)
  if [[ $lc == 0 ]]; then
    cat contigs2/${contig}/${contig}.fasta >> cur.v3.racon_polished2.fasta
  else
    cat contigs2/${contig}/racon_${contig}.no_trimming.fasta >> cur.v3.racon_polished2.fasta
  fi
done
hometools seqsize cur.v3.racon_polished2.fasta > cur.v3.racon_polished2.fasta.chrsize


# Map the HiFI reads to assembly, also generate dotplot between polished contigs and GenomeBiology assembly
cd mapping
minimap2 -ax asm5 --eqx -t 60 \
  /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/variantCorrected_AND_kmer_masked_d100-kmer80bp-8_4_Currot_on_Original_v1.0.fasta \
  ../cur.v3.racon_polished2.fasta  > out.sam
drawdotplot('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/cur/mapping/out.sam', out='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/cur/mapping/out.pdf', height=20) # use function from SyRI2

indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/'
minimap2 -ax map-pb -t 60 --secondary=no \
  ../cur.v3.racon_polished2.fasta \
  ${indir}/haplotype-cur.fasta.gz \
  ${indir}/haplotype-unknown.fasta.gz \
| samtools sort -@ 60 -O BAM - \
> cur.v3.racon_polished2.sorted.bam
samtools index -@60 cur.v3.racon_polished2.sorted.bam



# Racon third run polishing polished pilon contigs directly after removing alignments with mapping quality <10 and removing reads with clippings
# Align HiFi-reads for second RACON run
samtools view -@60 -h -q10 cur.v3.pilon_polished.sorted.bam \
| awk '$6 !~ /H|S/{print}' \
| samtools view -@60 -bS - \
> cur.v3.pilon_polished.q10.noclips.sorted.bam
samtools index -@60 cur.v3.pilon_polished.q10.noclips.sorted.bam

# Extract contig alignments and reads. NOTE: racon do not accept BAM and requires SAM as input
cat cur.v3.fasta.chrsize \
| cut -f1 \
| xargs -n 1 -P 60 -I {} bash -c '
    echo ${1}
    mkdir contigs3/${1}
    samtools view -h cur.v3.pilon_polished.q10.noclips.sorted.bam ${1}_pilon > contigs3/${1}/${1}.hifi.sorted.sam
    samtools fasta contigs3/${1}/${1}.hifi.sorted.sam | gzip > contigs3/${1}/${1}.hifi.fasta.gz
    gzip contigs3/${1}/${1}.hifi.sorted.sam
    hometools getchr --chrs ${1}_pilon -o contigs3/${1}/${1}.fasta cur.v3.pilon_polished.fasta
' -- {}


contigs=$(cat cur.v3.fasta.chrsize | cut -f1)
for contig in ${contigs[@]}; do
  cd contigs3/${contig}
  bsub -q short -n 3 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo racon.log -eo racon.err "
    racon -u --no-trimming -t 3 \
      ${contig}.hifi.fasta.gz \
      ${contig}.hifi.sorted.sam.gz \
      ${contig}.fasta > racon_${contig}.no_trimming.fasta
  "
  cd ../../
done

contigs=$(cat cur.v3.fasta.chrsize | cut -f1)
rm cur.v3.racon_polished3.fasta
for contig in ${contigs[@]}; do
  lc=$(wc -l contigs3/${contig}/racon_${contig}.no_trimming.fasta | cut -d' ' -f1)
  if [[ $lc == 0 ]]; then
    cat contigs3/${contig}/${contig}.fasta >> cur.v3.racon_polished3.fasta
  else
    cat contigs3/${contig}/racon_${contig}.no_trimming.fasta >> cur.v3.racon_polished3.fasta
  fi
done
hometools seqsize cur.v3.racon_polished3.fasta > cur.v3.racon_polished3.fasta.chrsize


# Map the HiFI reads to assembly, also generate dotplot between polished contigs and GenomeBiology assembly
cd mapping3
minimap2 -ax asm5 --eqx -t 60 \
  /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/variantCorrected_AND_kmer_masked_d100-kmer80bp-8_4_Currot_on_Original_v1.0.fasta \
  ../cur.v3.racon_polished3.fasta  > out.sam
drawdotplot('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/cur/mapping/out.sam', out='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/cur/mapping/out.pdf', height=20) # use function from SyRI2

indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/'
minimap2 -ax map-pb -t 60 --secondary=no \
  ../cur.v3.racon_polished3.fasta \
  ${indir}/haplotype-cur.fasta.gz \
  ${indir}/haplotype-unknown.fasta.gz \
| samtools sort -@ 60 -O BAM - \
> cur.v3.racon_polished3.sorted.bam
samtools index -@60 cur.v3.racon_polished3.sorted.bam


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


# Racon second run after removing alignments with mapping quality <10 and removing reads with clippings
# Align HiFi-reads for second RACON run
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/'
minimap2 -ax map-pb -t 60 --secondary=no \
  ora.v3.racon_polished.fasta \
  ${indir}/haplotype-ora.fasta.gz \
  ${indir}/haplotype-unknown.fasta.gz \
| samtools sort -@ 60 -O BAM - \
> ora.v3.racon_polished.sorted.bam
samtools view -h -q10 ora.v3.racon_polished.sorted.bam \
| awk '$6 !~ /H|S/{print}' \
| samtools view -bS - \
> ora.v3.racon_polished.q10.noclips.sorted.bam
samtools index -@60 ora.v3.racon_polished.q10.noclips.sorted.bam

# Extract contig alignments and reads. NOTE: racon do not accept BAM and requires SAM as input
cat ora.v3.fasta.chrsize \
| cut -f1 \
| xargs -n 1 -P 60 -I {} bash -c '
    echo ${1}
    mkdir contigs2/${1}
    samtools view -h ora.v3.racon_polished.q10.noclips.sorted.bam ${1}_pilon > contigs2/${1}/${1}.hifi.sorted.sam
    samtools fasta contigs2/${1}/${1}.hifi.sorted.sam | gzip > contigs2/${1}/${1}.hifi.fasta.gz
    gzip contigs2/${1}/${1}.hifi.sorted.sam
    hometools getchr --chrs ${1}_pilon -o contigs2/${1}/${1}.fasta ora.v3.racon_polished.fasta
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
  cd contigs2/${contig}
  bsub -q short -n 3 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo racon.log -eo racon.err "
    racon -u --no-trimming -t 3 \
      ${contig}.hifi.fasta.gz \
      ${contig}.hifi.sorted.sam.gz \
      ${contig}.fasta > racon_${contig}.no_trimming.fasta
  "
  cd ../../
done

contigs=$(cat ora.v3.fasta.chrsize | cut -f1)
rm ora.v3.racon_polished2.fasta
for contig in ${contigs[@]}; do
  lc=$(wc -l contigs2/${contig}/racon_${contig}.no_trimming.fasta | cut -d' ' -f1)
  if [[ $lc == 0 ]]; then
    cat contigs2/${contig}/${contig}.fasta >> ora.v3.racon_polished2.fasta
  else
    cat contigs2/${contig}/racon_${contig}.no_trimming.fasta >> ora.v3.racon_polished2.fasta
  fi
done
hometools seqsize ora.v3.racon_polished2.fasta > ora.v3.racon_polished2.fasta.chrsize


# Map the HiFI reads to assembly, also generate dotplot between polished contigs and GenomeBiology assembly
cd mapping
minimap2 -ax asm5 --eqx -t 60 \
  /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/variantCorrected_AND_kmer_masked_d100-kmer80bp-8_4_OrangeRed_on_Original_v1.0.fasta \
  ../ora.v3.racon_polished2.fasta  > out.sam
drawdotplot('out.sam', height=20) # use function from SyRI2

indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/'
minimap2 -ax map-pb -t 60 --secondary=no \
  ../ora.v3.racon_polished2.fasta \
  ${indir}/haplotype-ora.fasta.gz \
  ${indir}/haplotype-unknown.fasta.gz \
| samtools sort -@ 60 -O BAM - \
> ora.v3.racon_polished2.sorted.bam
samtools index -@60 ora.v3.racon_polished2.sorted.bam




# Racon third run polishing polished pilon contigs directly after removing alignments with mapping quality <10 and removing reads with clippings
# Align HiFi-reads for second RACON run
samtools view -@60 -h -q10 ora.v3.pilon_polished.sorted.bam \
| awk '$6 !~ /H|S/{print}' \
| samtools view -@60 -bS - \
> ora.v3.pilon_polished.q10.noclips.sorted.bam
samtools index -@60 ora.v3.pilon_polished.q10.noclips.sorted.bam

# Extract contig alignments and reads. NOTE: racon do not accept BAM and requires SAM as input
cat ora.v3.fasta.chrsize \
| cut -f1 \
| xargs -n 1 -P 60 -I {} bash -c '
    echo ${1}
    mkdir contigs3/${1}
    samtools view -h ora.v3.pilon_polished.q10.noclips.sorted.bam ${1}_pilon > contigs3/${1}/${1}.hifi.sorted.sam
    samtools fasta contigs3/${1}/${1}.hifi.sorted.sam | gzip > contigs3/${1}/${1}.hifi.fasta.gz
    gzip contigs3/${1}/${1}.hifi.sorted.sam
    hometools getchr --chrs ${1}_pilon -o contigs3/${1}/${1}.fasta ora.v3.pilon_polished.fasta
' -- {}


contigs=$(cat ora.v3.fasta.chrsize | cut -f1)
for contig in ${contigs[@]}; do
  cd contigs3/${contig}
  bsub -q short -n 3 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo racon.log -eo racon.err "
    racon -u --no-trimming -t 3 \
      ${contig}.hifi.fasta.gz \
      ${contig}.hifi.sorted.sam.gz \
      ${contig}.fasta > racon_${contig}.no_trimming.fasta
  "
  cd ../../
done

contigs=$(cat ora.v3.fasta.chrsize | cut -f1)
rm ora.v3.racon_polished3.fasta
for contig in ${contigs[@]}; do
  lc=$(wc -l contigs3/${contig}/racon_${contig}.no_trimming.fasta | cut -d' ' -f1)
  if [[ $lc == 0 ]]; then
    cat contigs3/${contig}/${contig}.fasta >> ora.v3.racon_polished3.fasta
  else
    cat contigs3/${contig}/racon_${contig}.no_trimming.fasta >> ora.v3.racon_polished3.fasta
  fi
done
hometools seqsize ora.v3.racon_polished3.fasta > ora.v3.racon_polished3.fasta.chrsize


# Map the HiFI reads to assembly, also generate dotplot between polished contigs and GenomeBiology assembly
cd mapping3
minimap2 -ax asm5 --eqx -t 60 \
  /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/variantCorrected_AND_kmer_masked_d100-kmer80bp-8_4_OrangeRed_on_Original_v1.0.fasta \
  ../ora.v3.racon_polished3.fasta  > out.sam
drawdotplot('/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/ora/mapping/out.sam', out='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/ora/mapping/out.pdf', height=20) # use function from SyRI2

indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/'
minimap2 -ax map-pb -t 60 --secondary=no \
  ../ora.v3.racon_polished3.fasta \
  ${indir}/haplotype-ora.fasta.gz \
  ${indir}/haplotype-unknown.fasta.gz \
| samtools sort -@ 60 -O BAM - \
> ora.v3.racon_polished3.sorted.bam
samtools index -@60 ora.v3.racon_polished3.sorted.bam


################################################################################
##################### Assembly assessment using Merqury ########################
################################################################################

# Create meryl database for the rojo-passion reads
## Merge Reads
cat *_A_*_R1_*fastq.gz > merged_A_R1.fastq.gz &
cat *_A_*_R2_*fastq.gz > merged_A_R2.fastq.gz &
cat *_B_*_R1_*fastq.gz > merged_B_R1.fastq.gz &
cat *_B_*_R2_*fastq.gz > merged_B_R2.fastq.gz &
cat *_C_*_R1_*fastq.gz > merged_C_R1.fastq.gz &
cat *_C_*_R2_*fastq.gz > merged_C_R2.fastq.gz &
cat *_D_*_R1_*fastq.gz > merged_D_R1.fastq.gz &
cat *_D_*_R2_*fastq.gz > merged_D_R2.fastq.gz &

cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/rojo_passion/'
adapter='/srv/netscratch/dep_mercier/grp_schneeberger/projects/hyper_co/data/reads/shqMergedAdapters_Primers_representative_rc.fa'
meryl='/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin_jose/meryl/build/bin/meryl'
for s in A B C D; do
  cd $cwd
  bsub -q multicore40 -n 40 -R "rusage[mem=50000] span[hosts=1]" -M 60000 -oo ${s}.log -eo ${s}.err "
    # Trim reads
    skewer -r 0.1 -d 0.05 -k 8 -q 20 -l 75 -m pe \
      -t 40 \
      -x $adapter \
      -z -o ${s}_ql \
      merged_${s}_R1.fastq.gz \
      merged_${s}_R2.fastq.gz
    $meryl k=21 count threads=40 memory=100 output ${s}_R1.meryl ${s}_ql-trimmed-pair1.fastq.gz &
    $meryl k=21 count threads=40 memory=100 output ${s}_R2.meryl ${s}_ql-trimmed-pair2.fastq.gz &
    "
done

# Combine Meryl databases
$meryl union-sum output rp.meryl A_R1.meryl A_R2.meryl B_R1.meryl B_R2.meryl C_R1.meryl C_R2.meryl D_R1.meryl D_R2.meryl

# Merqury assessment
export MERQURY=/srv/netscratch/dep_mercier/grp_schneeberger/software/merqury-1.1/

# Get parental hapmer database
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/merqury_assessment/parental_hapmers
cur_meryl=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/currot/currot.meryl
ora_meryl=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/orangered/orangered.meryl
rp_meryl=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/rojo_passion/rp.meryl
nohup $MERQURY/trio/hapmers.sh $cur_meryl $ora_meryl $rp_meryl &

# Analyse cur and ora assemblies
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/merqury_assessment/
cur_meryl=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/merqury_assessment/parental_hapmers/currot.inherited.gt18.meryl
ora_meryl=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/merqury_assessment/parental_hapmers/orangered.inherited.gt20.meryl
rp_meryl=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/rojo_passion/rp.meryl
cur_asm=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/cur/cur.v3.racon_polished.fasta
ora_asm=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/ora/ora.v3.racon_polished.fasta

cd $cwd
nohup $MERQURY/merqury.sh $rp_meryl $cur_meryl $ora_meryl $cur_asm $ora_asm rp &

# Cur phasing accuracy: 99.98%
cat rp.hapmers.count | grep 'cur.v3' | awk '{g+=$3; b+=$4} END {print g,b,b/(b+g),g/(b+g)}'
# Cur phasing accuracy: 99.99%
cat rp.hapmers.count | grep 'ora.v3' | awk '{g+=$4; b+=$3} END {print g,b,b/(b+g),g/(b+g)}'

# Analyse cur and ora assemblies RUN2
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/merqury_assessment/run2/
cur_meryl=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/merqury_assessment/parental_hapmers/currot.inherited.gt18.meryl
ora_meryl=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/merqury_assessment/parental_hapmers/orangered.inherited.gt20.meryl
rp_meryl=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/rojo_passion/rp.meryl
cur_asm=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/cur/cur.v3.racon_polished2.fasta
ora_asm=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/ora/ora.v3.racon_polished2.fasta

cd $cwd
nohup $MERQURY/merqury.sh $rp_meryl $cur_meryl $ora_meryl $cur_asm $ora_asm rp &

# Cur phasing accuracy: 99.98%
cat rp.hapmers.count | grep 'cur.v3' | awk '{g+=$3; b+=$4} END {print g,b,b/(b+g),g/(b+g)}'
# Cur phasing accuracy: 99.99%
cat rp.hapmers.count | grep 'ora.v3' | awk '{g+=$4; b+=$3} END {print g,b,b/(b+g),g/(b+g)}'


# Analyse cur and ora assemblies RUN3
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/merqury_assessment/run3/
cur_meryl=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/merqury_assessment/parental_hapmers/currot.inherited.gt18.meryl
ora_meryl=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/merqury_assessment/parental_hapmers/orangered.inherited.gt20.meryl
rp_meryl=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/rojo_passion/rp.meryl
cur_asm=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/cur/cur.v3.racon_polished3.fasta
ora_asm=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/ora/ora.v3.racon_polished3.fasta

cd $cwd
nohup $MERQURY/merqury.sh $rp_meryl $cur_meryl $ora_meryl $cur_asm $ora_asm rp &

# Cur phasing accuracy: 99.98%
cat rp.hapmers.count | grep 'cur.v3' | awk '{g+=$3; b+=$4} END {print g,b,b/(b+g),g/(b+g)}'
# Cur phasing accuracy: 99.99%
cat rp.hapmers.count | grep 'ora.v3' | awk '{g+=$4; b+=$3} END {print g,b,b/(b+g),g/(b+g)}'


################################################################################
######################### SALSA assembly correction  ###########################
################################################################################

# Split HiC fastq.gz file for allowing faster processing
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_HiC/
nohup zcat merged_R1.fastq.gz | split -l 10000000 --numeric-suffixes=10 - merged_R1_split &
nohup zcat merged_R2.fastq.gz | split -l 10000000 --numeric-suffixes=10 - merged_R2_split &
ls merged_R*_split* | xargs -n 1 -I {} -P 60 gzip {}

## Running SALSA on separate haplotypes to fix contig errors https://github.com/marbl/SALSA/issues/71
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/SALSA_contig_check/


# Get diploid genome and index
cd $cwd
cat ../../polishing/cur/cur.v3.racon_polished3.fasta \
| sed 's/^>u/>cur_u/g'  \
| sed '/^>/ s/ .*//' \
> cur/cur.v3.racon_polished.fasta
cat ../../polishing/ora/ora.v3.racon_polished3.fasta \
| sed 's/^>u/>ora_u/g' \
| sed '/^>/ s/ .*//' \
> ora/ora.v3.racon_polished.fasta

cp /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/cur/cur_unk.p_utg.noseq.gfa cur/
cp /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assemble_phased_reads/hifiasm_assembly/ora/ora_unk.p_utg.noseq.gfa ora/

# Map reads
# step I.1. align Hi-C reads: total: 99,858,814+91,335,467 = 191,194,281 read pairs
#                                    ~ 57.35828 Gb
#                                    ~ 238.9929x/240 Mb ~ 119.4964x/480 Mb = 119.4964x/haploid.genome

# step I.2. Filtering SAM file
#           shq note: according to Nan Wang, they made it with DpnII, similar to Mboi, so "MBOI" should be the one for selecting sequence GATC!
#       see: https://international.neb.com/faqs/0001/01/01/what-s-the-difference-between-dpni-dpnii-mboi-and-sau3ai
#           FAQ: What's the difference between DpnI, DpnII, MboI, and Sau3AI?
#              They all recognize the same sequence but have different methylation sensitivities.
#              DpnI will only cleave fully-adenomethylated dam sites and hemi-adenomethylated dam sites 60X more slowly.
#              DpnII and # MboI share methylation sensitivity and cleave dam sites which lack adenomethylation;
#                  each is blocked by complete dam methylation and probably by hemi-methylation as well.
#              Sau3AI will cleave all dam # sites regardless of adenomethylation
#                  but is completely blocked by methylated cytosines within the dam sequence.
#               see also: https://github.com/tangerzhang/ALLHiC/issues/52
#

# run salsa2
# set up bashrc for salsa2 python2.7 and required: export PATH=/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin/anaconda/install/envs/salsa2_env/bin/:$PATH

indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_HiC/
allhic=/srv/netscratch/dep_mercier/grp_schneeberger/bin/ALLHiC/
for s in cur ora; do
  cd $cwd/$s
  a=$( seq 10 86 )
  a_list=''; for i in ${a[@]}; do a_list+="${i} "; done
  bsub -q bigmem -n 40 -R "rusage[mem=100000]" -R "span[hosts=1]" -M 110000 -oo bwa_${s}.log -eo bwa_${s}.err "
#    a=\$( seq 10 86 )
#    bwa index -a bwtsw ${s}.v3.racon_polished.fasta
#    samtools faidx ${s}.v3.racon_polished.fasta
#    bwa aln -t 40 ${s}.v3.racon_polished.fasta ${indir}/merged_R1.fastq.gz > ${s}_R1.sai
#    bwa aln -t 40 ${s}.v3.racon_polished.fasta ${indir}/merged_R2.fastq.gz > ${s}_R2.sai
    echo $a_list | xargs -d ' ' -i -P 25 bash -c '
      bwa sampe ${s}.v3.racon_polished.fasta \
        -P \
        ${s}_R1.sai ${s}_R2.sai \
        ${indir}/merged_R1_split\${1}.gz \
        ${indir}/merged_R2_split\${1}.gz \
      | samtools sort -@ 40 -O BAM - \
      > ${s}_split\${1}.sorted.bam
    ' -- {}
    samtools merge -O BAM -@40 ${s}.sorted.bam ${s}_split*.sorted.bam
    ${allhic}/scripts/PreprocessSAMs.pl ${s}.sorted.bam ${s}.v3.racon_polished.fasta MBOI;
    ${allhic}/scripts/filterBAM_forHiC.pl ${s}.sorted.REduced.paired_only.bam ${s}_clean.sam;
    samtools view -bt ${s}.v3.racon_polished.fasta.fai ${s}_clean.sam > ${s}_clean.bam;
    bamToBed -i ${s}_clean.bam | sort -k 4 > ${s}_clean.bed;
    export PATH=/srv/netscratch/dep_mercier/grp_schneeberger/bin/bin/anaconda/install/envs/salsa2_env/bin/:$PATH
    python /netscratch/dep_mercier/grp_schneeberger/bin/SALSA/run_pipeline.py \
        -a ${s}.v3.racon_polished.fasta \
        -l ${s}.v3.racon_polished.fasta.fai \
        -g ${s}_unk.p_utg.noseq.gfa \
        -b ${s}_clean.bed \
        -m yes \
        -e GATC \
        -o scaffolds_utgs_i10_noM \
        -s 240000000 \
        -i 10 \
        &> salsa2_run_pipeline_i10.log
  "
done

indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/
# Get plots of the breakpoint regions using samplot
for s in cur ora; do
  cd $cwd/$s
  while read c p; do
    echo $c $p
    chr=$(echo $c | sed "s/${s}_//g")
    echo $chr
    start=$((p-10000))
    end=$((p+10000))
    samplot plot \
      -b $indir${s}/mapping3/${s}.v3.racon_polished3.sorted.F2048.bam \
      -o ${chr}_$p.pdf \
      -c $chr \
      -s $start \
      -e $end \
      -q 10 \
      --min_event_size 1
  done < scaffolds_utgs_i10_noM/input_breaks
  pdfunite utg000*pdf contigs_regions.pdf
  rm utg000*pdf
done

# Contigs to be broken were selected manually by checking the coverages in the contigs_regions.pdf file
# List of chimeric contigs (to be broken) is in break_contigs_list
# The position for fragmentation is in scaffolds_utgs_i10_noM/input_breaks
# The fragmented sequence is in scaffolds_utgs_i10_noM/assembly.cleaned.fasta

indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/polishing/
# Get cleaned assembly with broken chimeric contigs
for s in cur ora; do
  cd $cwd/$s
  sed -i "s/^u/${s}_u/g" break_contigs_list
  hometools getchr -F break_contigs_list -v -o ${s}.v3.racon_polished.contig_break.fasta ${s}.v3.racon_polished.fasta
  frag_contig_id=$(cat break_contigs_list | tr '\n' '|')
  frag_contig_id=${frag_contig_id:0:-1}
  select_chr=$(grep -E $frag_contig_id scaffolds_utgs_i10_noM/assembly.cleaned.fasta | tr -d '>' | tr '\n' ' ' )
  hometools getchr --chrs $select_chr -o TMP.break_contigs.fasta scaffolds_utgs_i10_noM/assembly.cleaned.fasta
  cat TMP.break_contigs.fasta >> ${s}.v3.racon_polished.contig_break.fasta
  rm TMP.break_contigs.fasta
done

# Generate .hic file from SALSA scaffolds
# You need to have a jar file for Juicertools in order to generate .hic file.
# It can be found here (https://github.com/aidenlab/juicer/wiki/Download).
#     Once you download this file on your machine, set the path to the jar file (including the jar file itself)
#         in the convert.sh script to JUICER_JAR variable.
#     Once that's done, simply run convert.sh script as convert.sh SALSA_OUT_DIR,
#         where SALSA_OUT_DIR is the directory where the SALSA scaffolds and all the intermediate output resides.
#         Please make sure you have latest version of GNU core utils in order to make use of --parallel option in GNU sort.
#     After this script finishes execution, the .hic file will be generated in the SALSA_OUT_DIR as salsa_scaffolds.hic.
#         This file can be loaded to Juicebox and visualized.

/netscratch/dep_mercier/grp_schneeberger/bin/SALSA/convert.sh ./scaffolds &>convert.log&

################################################################################
##################### Group contigs into linkage groups ########################
################################################################################

# Use reference-based assembly methods to get the broader grouping of contigs
# Contigs which cannot be grouped to any chromosome would be added in all groups
# Scaffolding then can be done on individual contig-groups

# Check ragtag
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/contig_grouping/

# Generate genome by joining the linkage groups
cd $cwd
cat scaffolded_final_manual_upd_map_group*.fa > genetic_map_genome.fasta

indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/SALSA_contig_check/
for s in cur ora; do
  cd $cwd/$s
  ragtag.py scaffold \
    -q 10 \
    -i 0.5 \
    -o genmap_grouping \
    -t 60 \
    ../genetic_map_genome.fasta \
    $indir/${s}/${s}.v3.racon_polished.contig_break.fasta
  lgs=( $(grep '>' ../genetic_map_genome.fasta | tr -d '>') )
  for lg in ${lgs[@]}; do
    lg_chrs=$(grep $lg genmap_grouping/ragtag.scaffolds.agp | grep $s | cut -f6 | tr '\n' ' ')
    hometools getchr --chrs $lg_chrs -o ${lg}_contigs.fasta $indir/${s}/${s}.v3.racon_polished.contig_break.fasta
  done
done


################################################################################
########################### Hi-C based scaffolding  ############################
################################################################################

cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/

# Get diploid genome and index
cat ../polishing/cur/cur.v3.racon_polished3.fasta |  sed 's/^>u/>cur_u/g' > rp_dip.fasta
cat ../polishing/ora/ora.v3.racon_polished3.fasta |  sed 's/^>u/>ora_u/g' >> rp_dip.fasta
sed -i '/^>/ s/ .*//' rp_dip.fasta

indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_HiC/
allhic=/srv/netscratch/dep_mercier/grp_schneeberger/software/ALLHiC/
a=$( seq 10 86 )
a_list=''; for i in ${a[@]}; do a_list+="${i} "; done
bsub -q bigmem -n 40 -R "rusage[mem=100000]" -R "span[hosts=1]" -M 110000 -oo bwa_rp.log -eo bwa_rp.err "
  bwa index -a bwtsw rp_dip.fasta
  samtools faidx rp_dip.fasta
  bwa aln -t 40 rp_dip.fasta ${indir}/merged_R1.fastq.gz  > R1.sai
  bwa aln -t 40 rp_dip.fasta ${indir}/merged_R2.fastq.gz  > R2.sai
  echo $a_list | xargs -d ' ' -i -P 20 bash -c '
    bwa sampe rp_dip.fasta \
      -P \
      R1.sai R2.sai \
      ${indir}/merged_R1_split\${1}.gz \
      ${indir}/merged_R2_split\${1}.gz \
    | samtools sort -@ 40 -O BAM - \
    > rp_split\${1}.sorted.bam
  ' -- {}
  samtools merge -O BAM -@40 rp.sorted.bam rp_split*.sorted.bam
  ${allhic}/scripts/PreprocessSAMs.pl rp.sorted.bam rp_dip.fasta MBOI;
  ${allhic}/scripts/filterBAM_forHiC.pl rp.sorted.REduced.paired_only.bam rp_clean.sam;
  samtools view -bt rp_dip.fasta.fai rp_clean.sam > rp_clean.bam;
  ALLHiC_partition -b rp_clean.bam -r rp_dip.fasta -e GATC -k 32 -m 25
"

# Hi-C based scaffoling did not result in very good results, so contigs were grouped using genetic-map
# to scaffold linkage group, first I map the Hi-C reads to the diploid genome assembly (chimeric broken)
# Then all reads mapping the contigs of a linkage group would be selected and processed to do HiC based
# scaffolding
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/lg_read_set/
cd $cwd
cat ../SALSA_contig_check/cur/cur.v3.racon_polished.contig_break.fasta > rp_diploid.fasta
cat ../SALSA_contig_check/ora/ora.v3.racon_polished.contig_break.fasta >> rp_diploid.fasta

indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_HiC/
bsub -q ioheavy -n 40 -R "rusage[mem=100000]" -R "span[hosts=1]" -M 110000 -oo bwa_rp.log -eo bwa_rp. "
  bwa mem \
    -t 40 \
    -T 75 \
    -a \
    rp_diploid.fasta \
    ${indir}/merged_R1.fastq.gz \
    ${indir}/merged_R2.fastq.gz \
  | samtools sort -@ 40 -O BAM - \
  > rp_diploid.hic.sorted.sam
"
# Separate HiC reads into different read sets
bsub -q normal -n 1 -R "rusage[mem=10000]" -R "span[hosts=1]" -M 15000 -oo read_sep.log -eo read_sep.err "
  syripy /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/separate_reads_to_lg.py
  gzip */*/*fq
"

#           shq note: according to Nan Wang, they made it with DpnII, similar to Mboi, so "GATC" should be the one!
#           In molecular biology, the DpnII restriction endonuclease family is a family of restriction endonucleases
#                                 which includes DpnII from Diplococcus pneumoniae.
#                                 These enzymes recognise the double-stranded DNA unmethylated sequence GATC
#                                 and cleave before G-1, where it encompasses the full length of the protein.[1]

lgdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/contig_grouping/
readsdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/lg_read_set/
allhic=/srv/netscratch/dep_mercier/grp_schneeberger/software/ALLHiC/
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/lg_hic_scaffold/
for s in cur ora; do
  for i in {1..8}; do
    cd $cwd
    mkdir $s
    cd $s
    mkdir lg${i}
    cd lg${i}
    bsub -q short -n 1 -R "rusage[mem=10000]" -R "span[hosts=1]" -M 12500 -oo bwa_${s}_${i}.log -eo bwa_${s}_${i}.err "
      ln -sf ${lgdir}/${s}/final_manual_upd_map_group${i}_contigs.fasta .
      ln -sf ${readsdir}/${s}lgs/lg${i}/R1.fq.gz .
      ln -sf ${readsdir}/${s}lgs/lg${i}/R2.fq.gz .
      bwa index -a bwtsw final_manual_upd_map_group${i}_contigs.fasta
      samtools faidx final_manual_upd_map_group${i}_contigs.fasta
      bwa aln -t 5 final_manual_upd_map_group${i}_contigs.fasta R1.fq.gz  > R1.sai
      bwa aln -t 5 final_manual_upd_map_group${i}_contigs.fasta R2.fq.gz  > R2.sai
      bwa sampe final_manual_upd_map_group${i}_contigs.fasta \
        -P \
        R1.sai R2.sai \
        R1.fq.gz R2.fq.gz \
      | samtools sort -@ 5 -O BAM - \
      > lg${i}.sorted.bam
      # Allhic make_bed_around_RE_site_modified_MG.pl script is bugged and do not print the last contig.
      ${allhic}/scripts/PreprocessSAMs_modified_MG.pl lg${i}.sorted.bam final_manual_upd_map_group${i}_contigs.fasta MBOI;
      ${allhic}/scripts/filterBAM_forHiC.pl lg${i}.sorted.REduced.paired_only.bam lg${i}_clean.sam;
      samtools view -bt final_manual_upd_map_group${i}_contigs.fasta.fai lg${i}_clean.sam > lg${i}_clean.bam;
      ALLHiC_partition -b lg${i}_clean.bam -r final_manual_upd_map_group${i}_contigs.fasta -e GATC -k 1 -m 25;
      allhic extract lg${i}_clean.bam final_manual_upd_map_group${i}_contigs.fasta --RE GATC;
      allhic optimize lg${i}_clean.counts_GATC.1g1.txt lg${i}_clean.clm
      ALLHiC_build final_manual_upd_map_group${i}_contigs.fasta
      hometools seqsize groups.asm.fasta > chrn.list
      ALLHiC_plot lg${i}_clean.bam groups.agp chrn.list 200k pdf
      ALLHiC_plot lg${i}_clean.bam groups.agp chrn.list 200k png
    "
    echo $s $i
  done
done


cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/lg_hic_scaffold/cur
cat */groups.asm.fasta > cur.scaffold.v1.fasta
cat */groups.agp > cur.scaffold.v1.agp
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/lg_hic_scaffold/ora
cat */groups.asm.fasta > ora.scaffold.v1.fasta
cat */groups.agp >ora.scaffold.v1.agp

gbdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/
ragtagdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/contig_grouping/
hicdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/lg_hic_scaffold/
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/assembly_comparison/
for s in cur ora; do
  cd ${cwd}/$s
  gbref=$(ls ${gbdir}/${s}*.filtered.fasta)
  echo $gbref
#  minimap2 -ax asm5 --eqx -t10 $gbref ${ragtagdir}/${s}/genmap_grouping/ragtag.scaffolds.fasta > ref_ragtag.sam &
#  minimap2 -ax asm5 --eqx -t10 $gbref ${hicdir}/${s}/${s}.scaffold.v1.fasta > ref_hic.sam &
#  minimap2 -ax asm5 --eqx -t10 ${ragtagdir}/${s}/genmap_grouping/ragtag.scaffolds.fasta ${hicdir}/${s}/${s}.scaffold.v1.fasta > ragtag_hic.sam &
  for outfmt in pdf png; do
    syripy /biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri2/syri2/src/py/drawsamplot.py \
      --qagp ragtag.scaffolds.agp \
      -o ref_ragtag.${outfmt} \
      ref_ragtag.sam &
    syripy /biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri2/syri2/src/py/drawsamplot.py \
      --qagp ${s}.scaffold.v1.agp \
      -o ref_hic.${outfmt} \
      ref_hic.sam &
    syripy /biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri2/syri2/src/py/drawsamplot.py \
      --ragp ragtag.scaffolds.agp \
      --qagp ${s}.scaffold.v1.agp \
      -o ragtag_hic.${outfmt} \
      ragtag_hic.sam &
  done
done


################################################################################
########################### Get organelle genome  ##############################
################################################################################
# Get ungrouped contigs
asdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/contig_grouping/
for s in cur ora; do
  cd $asdir/$s
  chrs=$(grep -v 'final' genmap_grouping/ragtag.scaffolds.agp | grep -v '#' |  sort -n -k 3 -r | cut -f1)
  hometools getchr --chrs $chrs -o ungrouped.fasta genmap_grouping/ragtag.scaffolds.fasta
done

indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/organelle_sequence/
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/plastid_dna/
for s in cur ora; do
  cd $cwd/$s
  refs=( $(ls -d ${indir}*/ | xargs -n 1 basename) )
  rm plastid_contigs.txt
  for ref in ${refs[@]}; do
#    minimap2 -x asm5 -c --eqx -t4 ${indir}/${ref}/plastid.fasta ${asdir}/${s}/ungrouped.fasta > ${s}_${ref}.paf &
    python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_plastid_contigs.py \
    ${s}_${ref}.paf \
    >> plastid_contigs.txt
  done
  cat plastid_contigs.txt | cut -f1 | sort -u > plastid_contigs_unique.txt
done


################################################################################
######################## Manual assembly correction  ###########################
################################################################################
# CUR
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/lg_hic_scaffold/
for s in cur ora; do
  cd ${cwd}/$s
  cat lg*/groups.agp > all_groups.agp
done
# agp files were then manually edited

contigdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/SALSA_contig_check/
agpdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/hic_scaffolding/lg_hic_scaffold/
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/manual_curated_assembly/
for s in cur ora; do
  cd ${cwd}/$s
  python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_fasta_from_agp.py \
    ${agpdir}${s}/all_groups_manual_edited.agp.txt \
    ${contigdir}/${s}/${s}.v3.racon_polished.contig_break.fasta \
    -o ${s}_manually_curated.v1.fasta
done

gbdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/
ragtagdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/contig_grouping/
for s in cur ora; do
  cd ${cwd}/$s
  gbref=$(ls ${gbdir}/${s}*.filtered.fasta)
  echo $gbref
  minimap2 -ax asm5 --eqx -t10 $gbref ${s}_manually_curated.v1.fasta > ref_man_cur_v1.sam
  minimap2 -ax asm5 --eqx -t10 ${ragtagdir}/${s}/genmap_grouping/ragtag.scaffolds.fasta ${s}_manually_curated.v1.fasta > ragtag_man_cur_v1.sam
  for outfmt in pdf png; do
    syripy /biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri2/syri2/src/py/drawsamplot.py \
      -o ref_man_cur_v1.${outfmt} \
      ref_man_cur_v1.sam &
    syripy /biodata/dep_mercier/grp_schneeberger/projects/SynSearch/scripts/python/syri2/syri2/src/py/drawsamplot.py \
      -o ragtag_man_cur_v1.${outfmt} \
      ragtag_man_cur_v1.sam &
  done
done

# Get unplaced and plastid contigs
plastiddir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/plastid_dna/
for s in cur ora; do
  cd ${cwd}/$s
  hometools getchr -F ${plastiddir}/${s}/plastid_contigs_unique.txt \
    -o plastid_contigs.fasta \
    ${ragtagdir}/${s}/ungrouped.fasta
  hometools getchr -F ${plastiddir}/${s}/plastid_contigs_unique.txt \
    -v \
    -o unplaced_contigs.fasta \
    ${ragtagdir}/${s}/ungrouped.fasta
done


################################################################################
###################### Genetic Mapp based scaffolding  #########################
################################################################################

# NOT USED IN FINAL ANALYSIS
# Using reference-based assembly methods

# Check ntjoin
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/geneticmap_scaffolding/ntjoin/
sed -i '/^>/ s/ .*//' cur.v3.racon_polished3.fasta
/srv/netscratch/dep_mercier/grp_schneeberger/private/manish/toolbox/Parser/FASTA/fasta_oneliner.pl currot.v1.1.fasta > currot.v1.1.oneline.fasta
/srv/netscratch/dep_mercier/grp_schneeberger/private/manish/toolbox/Parser/FASTA/fasta_oneliner.pl cur.v3.racon_polished3.fasta > cur.v3.racon_polished3.oneline.fasta
ntJoin assemble target=cur.v3.racon_polished3.oneline.fasta target_weight=1 references='currot.v1.1.oneline.fasta' reference_weights='2' k=32 w=500 t=5 agp=True no_cut=True
## ntJoin didnt result in good results

# Check ragtag
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hifi_assembly/assembly_polishing_scaffolding/geneticmap_scaffolding/ratag/
/srv/netscratch/dep_mercier/grp_schneeberger/software/RagTag-1.0.2/ragtag.py scaffold


