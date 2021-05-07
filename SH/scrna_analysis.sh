## Step 1: Make reference genome index
indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/
cd $indir
#gffread cur.pasa_out.gff -T > cur.pasa_out.gtf # Using the PASA out GFF for now, will need to udpate it when the final annotation is ready
gffread cur.pasa_out.sort.protein_coding.3utr.gff3 -T > cur.pasa_out.sort.protein_coding.3utr.gtf
# can edit the gff file manually to add UTRs

# Cellranger index
nohup cellranger mkref --genome=cur --fasta cur.genome.v1.fasta --genes=cur.pasa_out.sort.protein_coding.3utr.gtf --nthreads=20 &


# STAR Index
cd $indir
nohup STAR  --runThreadN 40 \
      --runMode genomeGenerate \
      --genomeSAindexNbases 12 \
      --genomeDir cur.genome.v1.fasta.star_index \
      --genomeFastaFiles cur.genome.v1.fasta \
      --sjdbGTFfile cur.pasa_out.gtf > np.star_index &


## Step 2: Perform barcode-correction and get cell count using cellranger count
indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scrna/bigdata/
refdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur/
cwd=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/get_cells/
cd $cwd
samples=( 'WT_1' 'WT_19' 'MUT_11_1' 'MUT_15' )
for sample in ${samples[@]}; do
#  cd ${indir}/${sample}
#  for f in $(ls *.fastq.gz); do
#    echo $f
#    mv $f  $(echo $f | sed -r 's/_H.{8}_/_/')
#  done
  cd $cwd
  mkdir $sample
  cd $sample
  sample_list=$(ls ${indir}/${sample}/*fastq.gz  | sed 's/^.*\///g' |sed 's/_S.*$//g' | sort -u | tr '\n' ',' | sed 's/,$//g')
  echo $sample_list

  bsub -q multicore20 -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo ${sample}.log -eo ${sample}.err "
    cellranger count --id=$sample \
      --fastqs=${indir}/${sample} \
      --transcriptome=$refdir \
      --sample=$sample_list \
      --expect-cells 3500 \
      --include-introns \
      --jobmode=lsf --maxjobs=10000 --jobinterval=1 --mempercore=10 \
      2>&1  > cellranger.log
  "

#  cellranger count --id=$sample \
#    --fastqs=${indir}/${sample} \
#    --transcriptome=$refdir \
#    --sample=$sample_list \
#    --expect-cells 3500 \
#    --localcores=15 \
#    --localmem=70 \
#    --include-introns \
#    2>&1  > cellranger.log \
#   &
  cd ..
done



cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/get_cells3/WT_1_reanalyse
cellranger reanalyze --id=WT_1_reanalyse \
  --matrix=../WT_1/WT_1/outs/filtered_feature_bc_matrix.h5 \
  --params=params.csv \
  --localcores=60 \
  --localmem=300



#############################################################
## Step 2: Align reads to genome index
#############################################################X
ref_idx='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta.star_index'
STAR --runThreadN 50 \
    --genomeDir $ref_idx \
    --readFilesIn /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/rojo_passion_rna/leaf_bulk/4031_D_run576_GCTAACTC_S32_L003_R1_001.fastq.gz,/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/rojo_passion_rna/leaf_bulk/4031_D_run576_GCTAACTC_S32_L004_R1_001.fastq.gz \
    --readFilesCommand zcat \
    --outStd BAM_SortedByCoordinate \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMprimaryFlag AllBestScore \
    --outFileNamePrefix rp_leaf \
    --quantMode GeneCounts \
    --twopassMode Basic \
    --alignIntronMax 20000 \
    --outFilterMismatchNmax 5 > rp_leaf.sorted.bam
