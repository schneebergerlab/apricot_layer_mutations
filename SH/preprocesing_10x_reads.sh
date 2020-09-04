## Step 1: Index current Currot and OrangeRed assemblies to be used as reference genomes for barcode correction

<<Comment
indir=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/

cd $indir
cellranger-dna mkref currot.v1.1.fasta currot_contig_defs.json
cellranger-dna mkref orangeRed.v1.1.fasta orangered_contig_defs.json


## Step 2: Use the indexed genomes and run the cellranger-dna CNV identification pipeline. One of the internediate steps for this pipeline is to barcode correction. We can use the final bam file to get the corrected barcodes.
## I do this step for both reference genomes (CUR and OR) and for all test libraries that includes:
## Two WT libraries: 4747_A and 4747_B
## Two MUT libraries: 4747_C and 4747_D
## First reads from different 10x barcodes (but corresponding to same libraries) are merged
## Then all reads smaller than 50bps were filtered out as cellranger-dna cannot work with smaller reads and there were naked barcode sequenced as reads. (Question: why did naked barcode ended up getting sequenced?)


run_barcode_correction () {
    bsub -q multicore40 -n 40 -R "span[hosts=1] rusage[mem=30000]" -M 30000 -oo $1.log -eo $1.err "
        cellranger-dna cnv --id=$1 --reference=$2 --fastq=$3 --sample=$4 --localcores=40 --localmem=30
        "
}        

refcur='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/refdata-currot.v1.1/'
refor='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/refdata-orangeRed.v1.1/'
wt='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/wt/'
mut='/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/mut/'


python=/netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/syri3.7/bin/python
filter=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/filter_short_reads.py

cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/wt/
#$python $filter 4747_A_merged_S1_L006_R1_001.fastq.gz 4747_A_merged_S1_L006_R2_001.fastq.gz 50
$python $filter 4747_B_merged_S5_L006_R1_001.fastq.gz 4747_B_merged_S5_L006_R2_001.fastq.gz 50

cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_scdna/mut/
$python $filter 4747_C_merged_S9_L006_R1_001.fastq.gz 4747_C_merged_S9_L006_R2_001.fastq.gz 50
$python $filter 4747_D_merged_S13_L006_R1_001.fastq.gz 4747_D_merged_S13_L006_R2_001.fastq.gz 50


cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/cur_wtA
run_barcode_correction cur_wtA $refcur $wt L50_4747_A_merged


cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/cur_wtB
run_barcode_correction cur_wtB $refcur $wt L50_4747_B_merged


cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/cur_mutC
run_barcode_correction cur_mutC $refcur $mut L50_4747_C_merged

cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/cur_mutD
run_barcode_correction cur_mutD $refcur $mut L50_4747_D_merged


cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/or_wtA
run_barcode_correction or_wtA $refor $wt L50_4747_A_merged


cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/or_wtB
run_barcode_correction or_wtB $refor $wt L50_4747_B_merged


cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/or_mutC
run_barcode_correction or_mutC $refor $mut L50_4747_C_merged


cd /srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/or_mutD
run_barcode_correction or_mutD $refor $mut L50_4747_D_merged
Comment

## Step 3: Separate reads from the corrected barcode bam file to separate folders corresponding to each barcode, get read count, and align them to the reference
T10xbam2fq=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/binaries/T10xbam2fq
asCellseparator=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/binaries/asCellseparator
curidx=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/initial_assembly_from_jose/currot.v1.1.fasta.mm2_Xsr.idx

get_cell_reads () {
    bsub -q normal -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo get_cell_reads.log -eo get_cell_reads.err "
            samtools sort -n $1 -o ${2}_${3}_RNsorted_bam.bam
            samtools view ${2}_${3}_RNsorted_bam.bam | $T10xbam2fq - ${2}_${3}
            mkdir barcodes        
            $asCellseparator 16 ${2}_${3}_fqfrom10xBam_bxCorrected_R1.fq.gz ${2}_${3}_fqfrom10xBam_bxCorrected_R2.fq.gz 10 2000000 barcodes/
            "
}

merge_fastqs () {
bsub -q normal -R "span[hosts=1] rusage[mem=10000]" -M 10000 -oo merge_fastqs.log -eo merge_fastqs.err " 
#    readlink -f ${1}/* > barcodes_list
#    while read r; do
#        cd \$r
#        bc=\$(basename \$r)
#        cat part*R1*.fq.gz > \${bc}_R1.fastq.gz
#        cat part*R2*.fq.gz > \${bc}_R2.fastq.gz
#    done < barcodes_list
    
    cd $2
    rm barcodes_read_count
    while read r; do
        bc=\$(basename \$r)
        cnt=\$(zgrep -c '@' \$r/\${bc}_R1.fastq.gz)
        echo \$bc \$cnt >> barcodes_read_count
    done < barcodes_list
    "
}

align_cells () {
    bsub -q normal -n 4 -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo align_cells.log -eo align_cells.err "
    xargs -a barcodes_list -n 1 -P 4 -I {} /bin/bash -c '
        cd {}
        bc=\$(basename {})
        minimap2 -ax sr \$2 \${bc}_R1.fastq.gz \${bc}_R2.fastq.gz | samtools view -O BAM - | samtools sort - > \${bc}.sorted.bam
        samtools index \${bc}.sorted.bam
        bamCoverage -b \${bc}.sorted.bam -of bedgraph -o \${bc}.bedgraph ' -- {} $1
"
}

#refs=( 'cur' 'or' )
refs=( 'cur' )
samples=( 'wtA' 'wtB' 'mutC' 'mutD' )

cwd=/srv/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/barcode_correction/

for r in ${refs[@]}; do
    for s in ${samples[@]}; do
        cd $cwd/${r}_${s}
    #    get_cell_reads ${r}_${s}/outs/possorted_bam.bam $r $s
        #merge_fastqs barcodes $cwd/${r}_${s}
        align_cells $curidx
    done
done



