## Testing whether we can observe long-range interactions between apricot genomic regions
## De novo mutations in the mutant samples that affect these long range interactions could be behind the phenotype

CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/hic_tad_calling/
INDIR=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_HiC/
CONFIG=/srv/netscratch/dep_mercier/grp_schneeberger/private/manish/configs/nfcore/hic/hic.config
cd $CWD
nextflow run nf-core/hic \
    -profile singularity \
    -r 1.3.0 \
    -c $CONFIG \
    --digestion 'dpnii' \
    --input '/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_HiC/*R{1,2}.fastq.gz' \
    --fasta /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta
