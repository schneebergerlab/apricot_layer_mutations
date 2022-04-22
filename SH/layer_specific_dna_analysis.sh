NOTE: FOR ALL ANALYSIS CURROT assembly would be used as the primary assembly

################################################################################
############# Pre-run: Get repeat regions in the currot assembly ###############
################################################################################
### Not required as done in the scdna_analysis.sh

################################################################################
############################ Pre-process reads #################################
################################################################################
hometools='/srv/netscratch/dep_mercier/grp_schneeberger/software/anaconda3_2021/envs/mgpy3.8/bin/python /srv/biodata/dep_mercier/grp_schneeberger/software/hometools/myUsefulFunctions.py'

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

# Copied the trimmed-read and alignments to the following folder and would now use that as CWD
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
    bsub -q ioheavy -n 8  -R "span[hosts=1] rusage[mem=5000]" -M 6000 -oo ${s}_bamrc.log -eo ${s}_bamrc.err  "
        $hometools pbamrc -n 8 -b 30 -q 10 -w 0 -S -I -f $refcur -l $CHRBED ${s}.deduped.bam bam_read_counts_b30_q10.bt2.txt

      # GET POSITIONS WITH AT LEAST THREE NON-REFERENCE BASES
      /netscratch/dep_mercier/grp_schneeberger/software/anaconda3/envs/syri3.8/bin/python /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/scripts/python/get_positions_with_low_ref_af.py bam_read_counts_b30_q10.bt2.txt
  "
done

grep 'MUT_11_1' /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mutations.regions > ${CWD}/mut_11_1.mutations.regions

MUTS=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/mutations.regions
MUTSNEW=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/high_cov_mutants_sorted.all_samples.unique.regions
for s in l1 l2 l3; do
    cd ${CWD}/mut_11_1_${s}

    ## Get read mapping depth histogram
#    cut -f 4 bam_read_counts_b30_q10.bt2.txt | sort -n | uniq -c | hometools plthist -o bam_read_counts_b30_q10.bt2.mapping_depth.hist.pdf -x Mapping Depth -y Frequency -t ${s}_bt2 -xlim 0 400 &

    # Get allele-frequency at the mutation positions in MUT_11_1
#    hometools pbamrc -b 30 -q 10 -w 0 -f $refcur -I -n 1 -l ../mut_11_1.mutations.regions ${s}.deduped.bam read_count_at_leaf_mut_pos.txt &

    # Get allele-frequency at the mutation positions in all branches
#    hometools pbamrc -b 30 -q 10 -w 0 -f $refcur -I -n 1 -l $MUTS ${s}.deduped.bam read_count_at_leaf_mut_pos.all_branches.txt &

#    Get allele-frequency at the updated mutation positions
    hometools pbamrc -b 30 -q 10 -w 0 -f $refcur -I -n 1 -l $MUTSNEW ${s}.deduped.bam read_count_at_leaf_mut_pos.all_branches.updated.txt &

done

## Plot allele frequency at MUT_11_1 positions
layer_specific_dna_analysis.py -> plot_snp_af()
## Plot allele frequency at all positions
layer_specific_dna_analysis.py -> plot_snp_af_all_branch()

####################################################################
############ De-Novo mutation identification
####################################################################
# Get read counts candidate layer_specific SMs
cd /netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples
sort -k1,1 -k2,2n -o layer_SM_candidates.txt layer_SM_candidates.txt
for sample in 'WT_1' 'WT_19' 'MUT_15' 'MUT_11_1' ; do
    hometools pbamrc -n 4 -b 0 -q 0 -w 0 -I -f $refcur -l layer_SM_candidates.txt ../${sample}/${sample}.sorted.bt2.bam ${sample}.sm_candidate.read_count.txt &
    echo $sample
done
for sample in 'wt7' 'wt18' 'mut4' 'mut11_2' ; do
    hometools pbamrc -n 4 -b 0 -q 0 -w 0 -I -f $refcur -l layer_SM_candidates.txt ../${sample}/${sample}.deduped.bam ${sample}.sm_candidate.read_count.txt &
    echo $sample
done


####################################################################
############ Gene-conversion identification
####################################################################
# Get read counts at SYRI snp positions
syrisnp=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/annotations/v1/haplodiff/syri_run/syri.snps.txt
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/
for s in l1 l2 l3 ; do
    cd ${CWD}/mut_11_1_${s}
    bsub -q multicore40 -n 40  -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ${s}_bamrc.log -eo ${s}_bamrc.err "
        $hometools pbamrc -b 30 -q 10 -w 0 -I -n 40 -f $refcur -l $syrisnp ${s}.deduped.bam ${s}.syri_snps.bamrc
    "
done




refcur=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/data/assemblies/hifi_assemblies/cur.genome.v1.fasta
CWD=/netscratch/dep_mercier/grp_schneeberger/projects/apricot_leaf/results/scdna/bigdata/variant_calling/layer_samples/
for sample in l1 l2 l3; do
    cd ${CWD}/mut_11_1_${sample}
    N=20
    for i in {1..8}; do
        snps="../../snps_split/snps_CUR${i}G.txt"
        bsub -q multicore40 -n $N -R "span[hosts=1] rusage[mem=10000]" -M 15000 -oo ${sample}_snps_pileup.log -eo ${sample}_snps_pileup.err "
            xargs -a $snps -P $N -I {} samtools mpileup -f $refcur -q 40 -E -Q 26 ${sample}.deduped.bam -r {} -O --output-QNAME > snps_CUR${i}G.pileup 2> garb
            sort -k2,2n -o snps_CUR${i}G.pileup snps_CUR${i}G.pileup
        "
    done
done




