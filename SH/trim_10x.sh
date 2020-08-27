trimmer='/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Plum/bin_Hequan_shared/T10X_barcode_trimmer'

## Each Step is one job. Other steps need to be commented out.

<<Comment
## Step 1:
### Trim library C reads
cwd='/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Plum/results/trimmed_10x/libC/'
indir='/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Plum/data/reads/10x/libC/'

cd $cwd
bsub -q normal -R "span[hosts=1] rusage[mem=1000]" -M 1000 -o C_trim.log -e C_trim.err "$trimmer ${indir}R1.fastq.gz ${indir}R2.fastq.gz";

### Trim library D reads
cwd='/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Plum/results/trimmed_10x/libD/'
indir='/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Plum/data/reads/10x/libD/'

cd $cwd
bsub -q normal -R "span[hosts=1] rusage[mem=1000]" -M 1000 -o d_trim.log -e D_trim.err "$trimmer ${indir}R1.fastq.gz ${indir}R2.fastq.gz";
Comment

## Step 2:

## Combine trimmed reads to get K-mers and size estimation
cwd='/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Plum/results/trimmed_10x/'
indir='/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Plum/results/trimmed_10x/'
cd $cwd

sample=plum
sizek=21

bsub -q multicore20 -n 20 -R "rusage[mem=30000]" -M 30000 -o log_jfcount_$sample -e err_jfcount_$sample "
    echo 1;
    zcat ${indir}/libC/trimmed_R1.fastq.gz ${indir}/libC/trimmed_R2.fastq.gz ${indir}/libD/trimmed_R1.fastq.gz ${indir}/libD/trimmed_R2.fastq.gz | jellyfish count /dev/fd/0  -C -o ${sample}_${sizek}mer_trimmed -m ${sizek} -t 20 -s 5G;
    jellyfish histo -h 300000 -o ${sample}_${sizek}mer_trimmed.histo ${sample}_${sizek}mer_trimmed;
"

## Step 3:
<<Comment
## Genome Size Estimation using findGSE
R 
library("findGSE")
findGSE(histo="plum_21mer_trimmed.histo", sizek=21, outdir=".", exp_hom=200)



## Step 4:
## Combine trimmed reads to one dataset
cwd='/srv/netscratch/dep_mercier/grp_schneeberger/projects/mutation_tree/Plum/results/trimmed_10x/'
cd $cwd
cat ${cwd}/libC/trimmed_R1.fastq.gz ${cwd}/libD/trimmed_R1.fastq.gz > plum_trimmed_R1.fastq.gz
cat ${cwd}/libC/trimmed_R2.fastq.gz ${cwd}/libD/trimmed_R2.fastq.gz > plum_trimmed_R2.fastq.gz
Comment






