## Step 1: Make reference genome index
<<Comment
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose
gffread Cur.v1.1.protein.coding.genes.gff -T > Cur.v1.1.protein.coding.genes.gtf
gffread Ora.v1.1.protein.coding.genes.gff -T > Ora.v1.1.protein.coding.genes.gtf

gffread Cur.v1.1.protein.coding.genes.manually_edited.gff -T > Cur.v1.1.protein.coding.genes.manually_edited.gtf

cellranger mkref --genome=cur --fasta currot.v1.1.fasta --genes=Cur.v1.1.protein.coding.genes.gtf --nthreads=20
cellranger mkref --genome=ora --fasta orangeRed.v1.1.fasta --genes=Ora.v1.1.protein.coding.genes.gtf --nthreads=20

cellranger mkref --genome=cur2 --fasta currot.v1.1.fasta --genes=Cur.v1.1.protein.coding.genes.manually_edited.gtf --nthreads=60 2>&1 > cur2_mkref.log

cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/v1.0
** Run the above commands again to generate the gtf file and then index the reference genomes **
Comment

## Step 2: Perform barcode-correction and get cell count using cellranger count
indir=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scrna/project_4801/
samples=( '4801_A_run672_SI-GA-C1' '4801_C_run672_SI-GA-E1' '4801_B_run672_SI-GA-D1'	'4801_D_run672_SI-GA-F1' )

<<Comment
## Run for v1.0 assembly
refdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/v1.0/cur/
cwd=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/v1.0/cur/
cd $cwd
for sample in ${samples[@]}; do
#  bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=30000]" -M 30000 -oo ${sample}.log -eo ${sample}.err "
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    cellranger count --id=$sample \
     --fastqs=$indir \
     --transcriptome=$refdir \
     --sample=$sample \
     --localcores=1 \
     --localmem=30 \
#  "
done
Comment


## Run for v1.1 assembly
refdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/cur/
cwd=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/v1.1/cur/
cd $cwd
for sample in ${samples[@]}; do
#  bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=30000]" -M 30000 -oo ${sample}.log -eo ${sample}.err "
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    cellranger count --id=$sample \
     --fastqs=$indir \
     --transcriptome=$refdir \
     --sample=$sample \
     --localcores=1 \
     --localmem=30 \
     &
#  "
done

## Run for RojoPassion diploid v1.0 assembly
refdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/v1.0/rp/
cwd=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/v1.0/rp/
cd $cwd
for sample in ${samples[@]}; do
#  bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=30000]" -M 30000 -oo ${sample}.log -eo ${sample}.err "
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    export JAVA_OPTIONS=-Xgcthreads 1
    cellranger count --id=$sample \
     --fastqs=$indir \
     --transcriptome=$refdir \
     --sample=$sample \
     --localcores=1 \
     --localmem=30 \
     &
#  "
done


############### Run for the second MUT_15 leaf ########
indir=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scrna/project_4917/
## Run for v1.1 assembly
refdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/cur/
cwd=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/v1.1/cur/
cd $cwd
samples=( '4917_A_run677_SI-GA-A2' )
for sample in ${samples[@]}; do
#  bsub -q multicore20 -n 1 -R "span[hosts=1] rusage[mem=30000]" -M 30000 -oo ${sample}.log -eo ${sample}.err "
    export MKL_NUM_THREADS=1
    export NUMEXPR_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    cellranger count --id=$sample \
     --fastqs=$indir \
     --transcriptome=$refdir \
     --sample=$sample \
     --localcores=50 \
     --localmem=100 \
     &
#  "
done

############### TEST Run for the bigdata using old CurV1.1 genome assembly ########
indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scrna/bigdata/
## Run for v1.1 assembly
refdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/cur/
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

  cellranger count --id=$sample \
    --fastqs=${indir}/${sample} \
    --transcriptome=$refdir \
    --sample=$sample_list \
    --localcores=12 \
    --localmem=50 \
    2>&1  > cellranger.log \
   &
  cd ..
done

############### TEST Run for the bigdata using old CurV1.1 genome assembly using Cellranger-5########
indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scrna/bigdata/
## Run for v1.1 assembly
refdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/cur/
cwd=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/get_cells2/
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

  cellranger count --id=$sample \
    --fastqs=${indir}/${sample} \
    --transcriptome=$refdir \
    --sample=$sample_list \
    --expect-cells 3500 \
    --localcores=12 \
    --localmem=50 \
    --include-introns \
    2>&1  > cellranger.log \
   &
  cd ..
done

############### TEST Run for the bigdata using old CurV1.1 genome assembly, Cellranger-5, and manually-edited annotations ########
indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scrna/bigdata/
## Run for v1.1 assembly
refdir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/cur2/
cwd=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/get_cells3/
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

  cellranger count --id=$sample \
    --fastqs=${indir}/${sample} \
    --transcriptome=$refdir \
    --sample=$sample_list \
    --localcores=15 \
    --localmem=70 \
    --include-introns \
    --expect-cells 3500 \
    2>&1  > cellranger.log \
   &
  cd ..
done
cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scrna/bigdata/get_cells3/WT_1_reanalyse
cellranger reanalyze --id=WT_1_reanalyse \
  --matrix=../WT_1/WT_1/outs/filtered_feature_bc_matrix.h5 \
  --params=params.csv \
  --localcores=60 \
  --localmem=300



# Align leaf bulk RNA to curV1.1
#############################################################
## Step 1: Index Reference Genome
#############################################################
STAR  --runThreadN 40 \
      --runMode genomeGenerate \
      --genomeDir /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/cur.v1.1_star_index \
      --genomeFastaFiles /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta \
      --sjdbGTFfile /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/Cur.v1.1.protein.coding.genes.manually_edited.gtf
#############################################################
## Step 2: Align reads to genome index
#############################################################X
ref_idx='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/cur.v1.1_star_index'
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


# Align leaf bulk RNA to curV1.1
#############################################################
## Step 1: Index Reference Genome
#############################################################
STAR  --runThreadN 40 \
      --runMode genomeGenerate \
      --genomeDir /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/cur.v1.1_old_anno_star_index \
      --genomeFastaFiles /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta \
      --sjdbGTFfile /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/Cur.v1.1.protein.coding.genes.gtf
#############################################################
## Step 2: Align reads to genome index
#############################################################X
ref_idx='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/cur.v1.1_old_anno_star_index'
STAR --runThreadN 50 \
    --genomeDir $ref_idx \
    --readFilesIn /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/rojo_passion_rna/leaf_bulk/4031_D_run576_GCTAACTC_S32_L003_R1_001.fastq.gz,/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/rojo_passion_rna/leaf_bulk/4031_D_run576_GCTAACTC_S32_L004_R1_001.fastq.gz \
    --readFilesCommand zcat \
    --outStd BAM_SortedByCoordinate \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMunmapped Within \
    --outSAMprimaryFlag AllBestScore \
    --outFileNamePrefix rp_leaf_oldanno \
    --quantMode GeneCounts \
    --twopassMode Basic \
    --alignIntronMax 20000 \
    --outFilterMismatchNmax 5 > rp_leaf_oldanno.sorted.bam
