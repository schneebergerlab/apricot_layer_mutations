## Step 1: Index current Currot and OrangeRed assemblies to be used as reference genomes for barcode correction

<<Comment
indir=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/

cd $indir
cellranger-dna mkref currot.v1.1.fasta currot_contig_defs.json
cellranger-dna mkref orangeRed.v1.1.fasta orangered_contig_defs.json
Comment

## Step 2: Use the indexed genomes and run the cellranger-dna CNV identification pipeline. One of the internediate steps for this pipeline is to barcode correction. We can use the final bam file to get the corrected barcodes.
## I do this step for both reference genomes (CUR and OR) and for all test libraries that includes:
##  Two WT libraries: 4747_A and 4747_B
##  Two MUT libraries: 4747_C and 4747_D

run_barcode_correction () {
    bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=30000]" -M 30000 -o $1.log -e $1.err "
        cellranger-dna cnv --id=$1 /
            --reference=$2 /
            --fastq=$3 /
            --sample=$4
            --localcores=40
            --localmem=30"
}        

refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/refdata-currot.v1.1/'
refor='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/refdata-orangeRed.v1.1/'
wt='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/wt/'
mut='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/mut/'

cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/cur_wtA
run_barcode_correction cur_wtA $refcur $wt 4747_A_merged


cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/cur_wtB
run_barcode_correction cur_wtB $refcur $wt 4747_B_merged


cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/cur_mutC
run_barcode_correction cur_mutC $refcur $mut 4747_C_merged

cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/cur_mutD
run_barcode_correction cur_mutD $refcur $mut 4747_D_merged


cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/or_wtA
run_barcode_correction or_wtA $refor $wt 4747_A_merged


cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/or_wtB
run_barcode_correction or_wtB $refor $wt 4747_B_merged


cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/or_mutC
run_barcode_correction or_mutC $refor $mut 4747_C_merged


cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/or_mutD
run_barcode_correction or_mutD $refor $mut 4747_D_merged





