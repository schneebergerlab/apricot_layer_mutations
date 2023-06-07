# Link data to working folder
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_isoseq/
cd $CWD
ln -s /srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4876/4876_A_run545_CCS.bam WT_1.iso_seq.ccs.bam
ln -s /srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4876/4876_B_run545_CCS.bam WT_19.iso_seq.ccs.bam
ln -s /srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4876/4876_C_run545_CCS.bam MUT_11_1.iso_seq.ccs.bam
ln -s /srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_4954/4954_A_run545_CCS.bam MUT_15.iso_seq.ccs.bam

################################################################################
############ Cluster Iso-seq reads into individual barcodes ####################
################################################################################

## Select reads having good primer orientation
### Lima from syri3.8 environment
SAMPLES=("MUT_11_1" "MUT_15" "WT_1" "WT_19")
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_isoseq/
# /srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/lima
# Library design -s leads to no good BAM identification
# Recommended workflow: https://isoseq.how/getting-started.html
wget https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-10x_barcodes/3M-february-2018-REVERSE-COMPLEMENTED.txt.gz

for sample in ${SAMPLES[@]}; do
    lima --isoseq --dump-clips --dump-removed -d ${sample}.iso_seq.ccs.bam 10x_primers.fasta ${sample}.iso_seq.bam &
done

declare -A sdict
sdict['WT_1']=85
sdict['WT_19']=85
sdict['MUT_11_1']=83
sdict['MUT_15']=80
# Select reads from barcodes and deduplicate reads (https://isoseq.how/umi/)
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    {
    ## Get BC and UMI tags for the reads
#    isoseq3 tag ${s}.iso_seq.5p--3p.bam ${s}.iso_seq.flt.bam --design T-12U-16B &
    ## Remove poly-A tails from the reads
#    isoseq3 refine --require-polya ${s}.iso_seq.flt.bam 10x_primers.fasta ${s}.iso_seq.flnc.bam &
    # Corect barcodes using the barcode white-list provided by 10x. Percentile cutoffs are selected based on the 'knee' in the plots. Consider increasing allowed mismatch between the RAW and corrected BC
#    isoseq3 correct --B 3M-february-2018-REVERSE-COMPLEMENTED.txt.gz -M 2 --method percentile --percentile ${sdict[$s]} ${s}.iso_seq.flnc.bam ${s}.iso_seq.flnc-bccorr.bam &
    # Get BC stats and check if the selected cutoff is correct
#    isoseq3 bcstats --method percentile --percentile ${sdict[$s]} -o ${s}.bcstats.tsv ${s}.iso_seq.flnc-bccorr.bam
#    python3 plot_knees.py -t ${s}.bcstats.tsv -o ${s}
    # Deduplicate reads based on UMIs
    samtools sort -@ 4 -t CB ${s}.iso_seq.flnc-bccorr.bam -o ${s}.iso_seq.flnc-bccorr.sorted.bam
#    isoseq3 groupdedup --keep-non-real-cells ${s}.iso_seq.flnc-bccorr.sorted.bam ${s}.iso_seq.flnc-bccorr.dedup.bam
    isoseq3 groupdedup ${s}.iso_seq.flnc-bccorr.sorted.bam ${s}.iso_seq.flnc-bccorr.dedup.bam
    } &
done
# The above pipeline resulted in the following number of reads:
# WT_1:  772070
# WT_19:  707564
# MUT_11_1:  534001
# MUT_15:  234369


## Get read count for each BC
indir='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_isoseq/'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/'
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $cwd
    mkdir -p $s; cd $s
    {
        samtools sort -t CB -O BAM ${indir}/${s}.iso_seq.flnc-bccorr.dedup.bam > ${s}.CB_sorted.bam
        samtools index ${s}.CB_sorted.bam
        samtools view ${s}.CB_sorted.bam  | grep -o -P 'CB:Z:[^\s]*' | uniq -c > cnt.txt
    } &
done

iso_seq_analysis.py --> get_iso_seq_stats()
# Separate Iso-seq reads to individual clusters
python.iso_seq_analysis.get_allele_freq_at_sm_pos_plot()
# Map iso-seq reads from individual clusters to cur reference and get read-count
# at SM positions
refcur='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
cwd='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/get_cells/'
for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    cd $cwd; cd $s
    clst=$(ls clstrs*reads.fa.gz)
    for c in ${clst[@]}; do
        cls=$(echo $c | sed 's/_reads.fa.gz//g')
        bsub -q multicore20 -n 10 -R "span[hosts=1] rusage[mem=20000]" -M 25000 -oo ${cls}.log -eo ${cls}.err "
            minimap2 -ax splice:hq -t 10 -R '@RG\tID:${cls}\tSM:1' --secondary=no -uf -C5 -O6,24 -B4 -Y $refcur $c \
            | samtools view -F 2048 -h - \
            | samtools sort -O BAM -@ 10 - \
            > ${cls}.iso_seq.bam
            samtools index -@ 10 ${cls}.iso_seq.bam
    "
    done
done

# Get readcount at all SM positions
## convert SM positions to regions files usable with bam-readcount
indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/
fin=${indir}/all_sm_in_all_samples.csv
tail +2 $fin \
| cut -f1,2 -d',' \
| awk -F ',' '{print $1"\t"$2"\t"$2}' > ${indir}/all_sm_in_all_samples.regions

# All sample somatic mutations
muts=${indir}/all_sm_in_all_samples.regions
for sample in WT_1 WT_19 MUT_15 MUT_11_1 ; do
    cd $cwd; cd $sample
    {
    for i in {1..12}; do
        hometools pbamrc -n 1 -b 0 -q 0 -w 0 -I -f $refcur -l $muts clstrs_${i}.iso_seq.bam clstrs_${i}.iso_seq.rc.txt
    done
    } &
done


# Full length transcript analysis. Trying to find allele-specifc expression (https://isoseq.how/classification/pigeon.html)
cwd=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/isoseq/allele_specific/
indir=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_isoseq/
refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta
curanno=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.pasa_out.sort.protein_coding.3utr.gtf
sqndir=/srv/netscratch/dep_mercier/grp_schneeberger/software/sqanti2_12_04_2023/SQANTI3-5.1.1/
echo /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/cur/RNAseq/4031.srt.bam > srbam.fofn
export PYTHONPATH=$PYTHONPATH:/srv/netscratch/dep_mercier/grp_schneeberger/software/cDNA_Cupcake/sequence/
export PYTHONPATH=$PYTHONPATH:/srv/netscratch/dep_mercier/grp_schneeberger/software/cDNA_Cupcake/

cd $cwd
# Sort and index the reference annotation files
#pigeon sort $curanno -o cur.annotation.sorted.gtf
#awk -i inplace -v FS=' ' '{print $0" gene_name "$12}' cur.annotation.sorted.gtf
#pigeon index cur.annotation.sorted.gtf

# Run sqanti_qc on the reference annotation file
python ${sqndir}/sqanti3_qc.py \
    cur.annotation.sorted.gtf \
    cur.annotation.sorted.gtf \
    $refcur \
    -o cur \
    -d . \
    --SR_bam srbam.fofn \
    -t 10 -n 4 --report skip --isoAnnotLite


for s in WT_1 WT_19 MUT_11_1 MUT_15; do
    {
        cd $cwd
        # Align reads found above the currot reference genome
#        pbmm2 align --preset ISOSEQ --sort -j 10 -J 4 ${indir}${s}.iso_seq.flnc-bccorr.dedup.bam $refcur ${s}.mapped.bam
        # Collapse redundant transcripts into unique isoforms based on exonic structures using isoseq collapse.
#        isoseq3 collapse ${s}.mapped.bam ${s}.collapsed.gff
        # Sort the transcript GFF file
#        pigeon sort ${s}.collapsed.gff -o ${s}.sorted.gff
        # Classify Isoforms (categories https://isoseq.how/classification/categories)
#        pigeon classify ${s}.sorted.gff cur.annotation.sorted.gtf $refcur --fl ${s}.collapsed.abundance.txt
        # Filter isoforms from the classification output.
#        pigeon filter ${s}_classification.txt --isoforms ${s}.sorted.gff -c 5 --min-cov 5

        # The classify and filter sub-commands from pigeon were resulting in quite many isoforms.
        # So testing the classification and filtering by SQANTI3

        # conda activate SQANTI3.env

        # Run qc
#        mkdir -p $s
        cd $s;
#        rm -r splits
#        tail +3 ../${s}.collapsed.gff > ${s}.collapsed.gtf
#        python ${sqndir}/sqanti3_qc.py \
#            ${s}.collapsed.gtf \
#            ../cur.annotation.sorted.gtf \
#            $refcur \
#            -o ${s} \
#            -d . \
#            -fl ../${s}.collapsed.abundance.txt \
#            --SR_bam ../srbam.fofn \
#            -t 10 -n 4 --report both --isoAnnotLite
        # Filter using default rules
#        python ${sqndir}/sqanti3_filter.py rules \
#            --isoAnnotGFF3 ${s}.gff3 \
#            --isoforms ${s}_corrected.fasta \
#            --gtf ${s}_corrected.gtf \
#            --faa ${s}_corrected.faa \
#            -o ${s} \
#            -d . \
#            ${s}_classification.txt
        # Filter using ML method
#        python ${sqndir}/sqanti3_filter.py ML \
#            --isoAnnotGFF3 ${s}.gff3 \
#            --isoforms ${s}_corrected.fasta \
#            --gtf ${s}_corrected.gtf \
#            --faa ${s}_corrected.faa \
#            -o ${s}.ML \
#            -d . \
#            ${s}_classification.txt
        # Rescue transcripts using sqanti3_rescue
        python ${sqndir}/sqanti3_rescue.py ml \
            --isoforms ${s}_corrected.fasta \
            --gtf ${s}.ML.filtered.gtf \
            --refGTF ../cur.annotation.sorted.gtf \
            -f $refcur \
            --refClassif ../cur_classification.txt \
            -o ${s}.rescue \
            -d . \
            -r randomforest.RData \
            ${s}.ML_MLresult_classification.txt
    } &
done
# The above pipeline resulted in following number of transcripts:
#MUT_11_1: 45390
#MUT_15: 22877
#WT_19: 57937
#WT_1: 61485

# Now check if there are any genes with clear difference in the transcriptome
iso_seq_analysis.py -> get_transcriptome_variants()


############ OLD STUFF BEFORE REANALYSIS OF ISO-SEQ READS ######################

/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/iso_seq_analysis.py

