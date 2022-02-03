Analyse the
NOTE: FOR ALL ANALYSIS CURROT assembly would be used as the primary assembly

################################################################################
############# Pre-run: Get repeat regions in the currot assembly ###############
################################################################################
### Not required as done in the scdna_analysis.sh

################################################################################
############################ Pre-process reads #################################
################################################################################
hometools='/netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/syri3.8/bin/python /srv/biodata/dep_mercier/grp_schneeberger/software/hometools/myUsefulFunctions.py'

### Get sequencing quality plots
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/rp_layer_seq/
INDIR=/srv/biodata/dep_mercier/grp_schneeberger/reads/Apricot/layer_specific/project_5363/X204SC21122744-Z01-F001/raw_data/
declare -A sample_libid
sample_libid['A']='mut_11_1_l1'
sample_libid['B']='mut_11_1_l2'
sample_libid['C']='mut_11_1_l3'

cd $CWD
for lib in A B C; do
    {
        cd $CWD
        mkdir ${sample_libid[$lib]}; cd ${sample_libid[$lib]}
        ln -s ${INDIR}/*_${lib}/*fq.gz .
        fastqc -t 10 *fq.gz
    } &
done

### Merge, trim and align reads to CUR genomes
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/reads/leaf_illumina/rp_layer_seq/
cd $CWD
#### Part below is not updated to the new folder structure
curidxbt2='/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta'
for s in l1 l2 l3; do
        cd ${CWD}/mut_11_1_${s}
        bsub -q multicore20 -n 20 -R "span[hosts=1] rusage[mem=20000]" -M 25000 -oo ${s}.log -eo ${s}.err "
            cat A5363*_1.fq.gz > merged_read_1.fq.gz
            cat A5363*_2.fq.gz > merged_read_2.fq.gz
            skewer -r 0.1 -d 0.05 -k 8 -q 20 -l 75 -m pe -t 15 -x \
            /srv/netscratch/dep_mercier/grp_schneeberger/projects/hyper_co/data/reads/shqMergedAdapters_Primers_representative_rc.fa \
            merged_read_1.fq.gz \
            merged_read_2.fq.gz \
            -z -o ${s}_ql

            bowtie2 --end-to-end \
                --very-sensitive \
                --threads 20 \
                -x $curidxbt2 \
                -1 ${s}_ql-trimmed-pair1.fastq.gz \
                -2 ${s}_ql-trimmed-pair2.fastq.gz \
                --rg-id ${s}\'\\t\'PL:ILLUMINA\'\\t\'SM:1 \
            | samtools view -b - \
            > ${s}.bam

            ## Mark duplicated reads
            samtools sort -@ 10 -n -O BAM ${s}.bam \
            | samtools fixmate -@ 10 -c -m -O BAM - - \
            | samtools sort -@ 10 -O BAM - \
            | samtools markdup -@ 10 -S -s -O BAM - - \
            | samtools sort -@ 10 -O BAM - > ${s}.DUPmarked.bam
            samtools index ${s}.DUPmarked.bam

            samtools view -G 1024 -O BAM ${s}.DUPmarked.bam > ${s}.deduped.bam
            samtools index ${s}.deduped.bam
        "
done

CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/
for s in l1 l2 l3; do
    cd ${CWD}/mut_11_1_${s}
    bsub -q normal -n 1 -R "span[hosts=1] rusage[mem=15000]" -M 20000 -oo ${s}_cov.log -eo ${s}_cov.err "
        $hometools bamcov -g ${s}.deduped.bam ${s}_vs_cur.deduped.cov
    "
done


####################################################################
############ Get read-counts at variant positions
####################################################################
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/
CHRBED=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.chrlen.bed
refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta
for s in l1 l2 l3; do
    cd ${CWD}/mut_11_1_${s}
    bsub -q normal -n 4  -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${s}_bamrc.log -eo ${s}_bamrc.err "
#        hometools pbamrc -n 4 -b 30 -q 10 -w 0 -S -f $refcur -l $CHRBED ${s}.deduped.bam bam_read_counts_b30_q10.bt2.txt

      # GET POSITIONS WITH AT LEAST THREE NON-REFERENCE BASES
      /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/syri3.8/bin/python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py bam_read_counts_b30_q10.bt2.txt
  "
done



## Get read mapping depth histogram
for s in l1 l2 l3; do
    cd ${CWD}/mut_11_1_${s}
    cut -d' ' -f 4 bam_read_counts_b30_q10.bt2.txt | sort -n | uniq -c | hometools plthist -o bam_read_counts_b30_q10.bt2.mapping_depth.hist.pdf -x Mapping Depth -y Frequency -t ${s}_bt2 -xlim 0 400 &
done

####################################################################
############ Gene-conversion identification
####################################################################
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/
HETPOSBED=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mitotic_recomb/strict_syn_snp_allele_readcount.depth450-650.af0.4-0.6.bed
refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta
for s in l1 l2 l3; do
    cd ${CWD}/mut_11_1_${s}
    bsub -q normal -n 4  -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${s}_bamrc_hetpos.log -eo ${s}_bamrc_hetpos.err "
        hometools pbamrc -n 4 -b 30 -q 10 -w 0 -f $refcur -l $HETPOSBED ${s}.deduped.bam hetpos.read_count.txt
  "
done

