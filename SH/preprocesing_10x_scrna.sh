## Step 1: Make reference genome index
<<Comment
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose
gffread Cur.v1.1.protein.coding.genes.gff -T > Cur.v1.1.protein.coding.genes.gtf
gffread Ora.v1.1.protein.coding.genes.gff -T > Ora.v1.1.protein.coding.genes.gtf

cellranger mkref --genome=cur --fasta currot.v1.1.fasta --genes=Cur.v1.1.protein.coding.genes.gtf --nthreads=20
cellranger mkref --genome=ora --fasta orangeRed.v1.1.fasta --genes=Ora.v1.1.protein.coding.genes.gtf --nthreads=20

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
